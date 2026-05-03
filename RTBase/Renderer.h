#pragma once

#include "Core.h"
#include "Sampling.h"
#include "Geometry.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"
#include "Scene.h"
#include "GamesEngineeringBase.h"
#include <thread>
#include <functional>
#include <vector>
#include <mutex>
#include <atomic>
#include <OpenImageDenoise/oidn.hpp>

class RayTracer
{
public:
	Scene* scene;
	GamesEngineeringBase::Window* canvas;
	Film* film;
	MTRandom *samplers;
	std::thread **threads;
	int numProcs;
	//10-1 准备（多线程分块-讲义没有）
	std::mutex filmSplatMutex;
	// 9-1 准备
	// AOV buffers: width * height * 3 (RGB float)
	std::vector<float> albedoBuffer;
	std::vector<float> normalBuffer;

	//9-2 准备
	// OIDN input: linear HDR color buffer (pre-tonemap)
	std::vector<float> colorBuffer;

	void init(Scene* _scene, GamesEngineeringBase::Window* _canvas)
	{
		scene = _scene;
		canvas = _canvas;
		film = new Film();
		film->init((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, new GaussianFilter());//选择滤波器

		SYSTEM_INFO sysInfo;
		GetSystemInfo(&sysInfo);
		numProcs = sysInfo.dwNumberOfProcessors;

		threads = new std::thread*[numProcs];
		samplers = new MTRandom[numProcs];

		// 4-4
		// Ensure each worker/thread gets a different RNG seed.
		// 为每个线程生成不同的随机种子
		std::random_device rd;
		const unsigned int baseSeed = rd();
		for (int i = 0; i < numProcs; ++i)
		{
			// 哈希偏移保证每个线程的种子互不相同
			const unsigned int seed = baseSeed ^ (0x9e3779b9u + (unsigned int)i * 0x85ebca6bu);
			samplers[i] = MTRandom(seed);
		}

		clear();
	}
	void clear()
	{
		film->clear();
	}
	//11-5原来只从环境贴图方向采一次样，现在同时从 BSDF 再采一次样，两者加权合并，各取所长
	// 5-2 当你随机挑出了一盏灯后，这一步要计算这盏灯是不是真的能照亮你正在看的这个物体表面。你需要检查光线有没有被别的物体挡住（可见性/阴影），并且要算上光线照射角度带来的能量衰减（也就是讲义中提到的几何项 / Geometry Term）
	Colour computeDirect(ShadingData shadingData, Sampler* sampler)
	{
		// 纯镜面材质不走这条“随机连线”的直接光照估计
		if (shadingData.bsdf->isPureSpecular() == true)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}

		float lightPMF = 0.0f;
		Light* light = scene->sampleLight(sampler, lightPMF);
		if (light == NULL || lightPMF <= 0.0f)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}

		Colour Le;
		float lightPDF = 0.0f;
		Vec3 lightSample = light->sample(shadingData, sampler, Le, lightPDF);
		if (lightPDF <= 0.0f)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}

		// 面光源：sample 返回的是光源表面点（面积域）
		if (light->isArea())
		{
			Vec3 toLight = lightSample - shadingData.x;
			float distSq = toLight.lengthSq();
			if (distSq <= (EPSILON * EPSILON))
			{
				return Colour(0.0f, 0.0f, 0.0f);
			}

			float dist = sqrtf(distSq);
			Vec3 wi = toLight / dist;

			float cosSurface = std::max(0.0f, Dot(shadingData.sNormal, wi));
			if (cosSurface <= 0.0f)
			{
				return Colour(0.0f, 0.0f, 0.0f);
			}

			Vec3 nLight = light->normal(shadingData, wi);
			float cosLight = std::max(0.0f, Dot(nLight, -wi));
			if (cosLight <= 0.0f)
			{
				return Colour(0.0f, 0.0f, 0.0f);
			}

			if (!scene->visible(shadingData.x, lightSample))
			{
				return Colour(0.0f, 0.0f, 0.0f);
			}

			Colour f = shadingData.bsdf->evaluate(shadingData, wi);
			float weight = (cosSurface * cosLight) / (distSq * lightPDF * lightPMF);
			return f * Le * weight;
		}

		// 非面光源（环境光）：MIS = Light Sampling + BSDF Sampling
		Colour Ld(0.0f, 0.0f, 0.0f);

		const Vec3 centre = use<SceneBounds>().sceneCentre;
		const float radius = use<SceneBounds>().sceneRadius;
		const float distToCentre = (shadingData.x - centre).length();
		const float tMax = distToCentre + radius + 1.0f;

		// [1] Light Sampling
		{
			Vec3 wiL = lightSample.normalize();
			const float pL = lightPDF;

			float cosTheta = std::max(0.0f, Dot(shadingData.sNormal, wiL));
			if (pL > 0.0f && cosTheta > 0.0f)
			{
				Vec3 farPoint = shadingData.x + (wiL * tMax);
				if (scene->visible(shadingData.x, farPoint))
				{
					float pB = shadingData.bsdf->PDF(shadingData, wiL);
					float denom = pL + pB;
					float wL = (denom > 0.0f) ? (pL / denom) : 0.0f;

					Colour f = shadingData.bsdf->evaluate(shadingData, wiL);
					Ld = Ld + (f * Le * (wL * cosTheta / (pL * lightPMF)));
				}
			}
		}

		// [2] BSDF Sampling
		{
			Colour fSample;
			float pB = 0.0f;
			Vec3 wiB = shadingData.bsdf->sample(shadingData, sampler, fSample, pB);

			float cosTheta = std::max(0.0f, Dot(shadingData.sNormal, wiB));
			if (pB > 0.0f && cosTheta > 0.0f)
			{
				Vec3 farPoint = shadingData.x + (wiB * tMax);
				if (scene->visible(shadingData.x, farPoint))
				{
					float pL = light->PDF(shadingData, wiB);
					float denom = pL + pB;
					float wB = (denom > 0.0f) ? (pB / denom) : 0.0f;

					Colour LeB = light->evaluate(wiB);
					Colour f = shadingData.bsdf->evaluate(shadingData, wiB);
					Ld = Ld + (f * LeB * (wB * cosTheta / (pB * lightPMF)));
				}
			}
		}

		return Ld;
	}
	// 5-3 让光线打中物体后，不仅算一下头顶上的灯光（上一步的直接光照），还要随机找个方向继续弹射出去，去收集“别的物体反射过来的光”。光线这样一路弹射形成的轨迹，就叫“路径（Path）”
	// 7-适配把 pathTrace 的“间接弹射”从均匀半球采样改为 BSDF::sample(...)，并对纯镜面分支做了特殊处理。
	Colour pathTrace(Ray& r, Colour& pathThroughput, int depth, Sampler* sampler, bool fromSpecular = true)
	{
		// Safety guard for recursion
		if (depth > 64)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}

		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);

		// Miss: accumulate environment/background
		if (shadingData.t >= FLT_MAX)
		{
			if (fromSpecular)
			{
				return pathThroughput * scene->background->evaluate(r.dir);
			}
			return Colour(0.0f, 0.0f, 0.0f);
		}

		Colour L(0.0f, 0.0f, 0.0f);

		// Hit emissive surface (avoid double counting with direct-light sampling)
		if (shadingData.bsdf->isLight())
		{
			if (fromSpecular)
			{
				L = L + (pathThroughput * shadingData.bsdf->emit(shadingData, shadingData.wo));
			}
		}

		// Add direct lighting at current bounce
		L = L + (pathThroughput * computeDirect(shadingData, sampler));

		// Indirect bounce: BSDF-driven sampling
		Colour f;
		float pdf = 0.0f;
		Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, f, pdf);

		// For invalid sample / impossible event
		if (pdf <= 0.0f)
		{
			return L;
		}

		Colour newThroughput;
		if (shadingData.bsdf->isPureSpecular())
		{
			// Mirror/Glass sample() already handles its own cosine-related term
			newThroughput = pathThroughput * f;
		}
		else
		{
			float cosTheta = std::max(0.0f, Dot(shadingData.sNormal, wi));
			if (cosTheta <= 0.0f)
			{
				return L;
			}
			newThroughput = pathThroughput * f * (cosTheta / pdf);
		}

		// Russian roulette termination
		const int rrStartDepth = 3;
		if (depth >= rrStartDepth)
		{
			float continueProb = std::max(0.05f, std::min(0.95f, newThroughput.Lum()));
			if (sampler->next() > continueProb)
			{
				return L;
			}
			newThroughput = newThroughput / continueProb;
		}

		// Clamp path throughput to suppress energy explosion / fireflies
		newThroughput.r = std::min(newThroughput.r, 10.0f);
		newThroughput.g = std::min(newThroughput.g, 10.0f);
		newThroughput.b = std::min(newThroughput.b, 10.0f);

		Ray nextRay;
		nextRay.init(shadingData.x + (wi * EPSILON), wi);

		bool nextFromSpecular = shadingData.bsdf->isPureSpecular();
		L = L + pathTrace(nextRay, newThroughput, depth + 1, sampler, nextFromSpecular);
		return L;
	}
	Colour direct(Ray& r, Sampler* sampler)
	{
		// Compute direct lighting for an image sampler here
		return Colour(0.0f, 0.0f, 0.0f);
	}
	Colour albedo(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return shadingData.bsdf->evaluate(shadingData, Vec3(0, 1, 0));
		}
		return scene->background->evaluate(r.dir);
	}
	Colour viewNormals(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		if (intersection.t < FLT_MAX)
		{
			ShadingData shadingData = scene->calculateShadingData(intersection, r);
			return Colour(fabsf(shadingData.sNormal.x), fabsf(shadingData.sNormal.y), fabsf(shadingData.sNormal.z));
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}
	//9-1（SpeedingItUp）渲染时额外收集 AOV 缓冲（法线 + 反照率）
	//8-end
	// 6-1 前我们用蒙特卡洛“扔骰子”法时，屏幕上每个像素都扔了同样次数的骰子（比如每个像素发 100 条光线）。但有些地方（比如平滑的白墙）很容易就算准了，有些地方（比如复杂的玻璃阴影）算很久还是有很多噪点（马赛克）。自适应采样就是教电脑“把好钢用在刀刃上”：聪明的分配计算资源，哪里噪点多，就往哪里多发光线
	// 2-4 把你上面写好的“虚拟胶片”、“色调映射”等零件，真正组装到你的主程序里，让它们跑起来
	void render()
	{
		// Keep film SPP semantics unchanged: one averaged contribution per pixel per frame.
		film->incrementSPP();

		const unsigned int width = film->width;
		const unsigned int height = film->height;
		const unsigned int pixelCount = width * height;

		//9-1
		// AOV buffer resize/reset (float RGB)
		albedoBuffer.assign(pixelCount * 3, 0.0f);
		normalBuffer.assign(pixelCount * 3, 0.0f);

		//9-2
		// OIDN input: linear HDR color (pre-tonemap)
		colorBuffer.assign(pixelCount * 3, 0.0f);
		std::vector<float> outputBuffer(pixelCount * 3, 0.0f);

		// 8-5 在现有路径追踪循环外，额外执行 light tracing，并按光路数归一化
		const int numLightPaths = (int)(width * height);
		const float lightPathScale = (numLightPaths > 0) ? (1.0f / (float)numLightPaths) : 0.0f;
		for (int i = 0; i < numLightPaths; ++i)
		{
			lightTrace(&samplers[0], lightPathScale);
		}

		const int blockSize = 16;
		const int pilotSPP = 2;
		const int maxSPPPerPixelThisFrame = 6;

		const int blocksX = (int)((width + blockSize - 1) / blockSize);
		const int blocksY = (int)((height + blockSize - 1) / blockSize);
		const int blockCount = blocksX * blocksY;

		std::vector<Colour> frameBuffer(width * height, Colour(0.0f, 0.0f, 0.0f));
		std::vector<float> blockSumLum(blockCount, 0.0f);
		std::vector<float> blockSumLumSq(blockCount, 0.0f);
		std::vector<int> blockPixelCount(blockCount, 0);
		std::vector<float> blockVariance(blockCount, 0.0f);
		std::vector<int> blockSPP(blockCount, pilotSPP);

		auto blockIndexOf = [blocksX, blockSize](int x, int y) -> int
			{
				const int bx = x / blockSize;
				const int by = y / blockSize;
				return by * blocksX + bx;
			};
		// 10-2 把现在两个嵌套 for 循环（遍历 y、x）的逻辑，改成按 tile 分发给多个线程执行。
		// 1) Pilot pass
		{
			const int workerCount = std::max(1, numProcs);
			std::atomic<int> nextTile(0);
			std::vector<std::thread> workers;
			workers.reserve((size_t)workerCount);

			for (int threadIndex = 0; threadIndex < workerCount; ++threadIndex)
			{
				workers.emplace_back([&, threadIndex]()
					{
						Sampler* sampler = &samplers[threadIndex];

						while (true)
						{
							const int tileId = nextTile.fetch_add(1);
							if (tileId >= blockCount)
							{
								break;
							}

							const int tileX = tileId % blocksX;
							const int tileY = tileId / blocksX;

							const unsigned int xBegin = (unsigned int)(tileX * blockSize);
							const unsigned int yBegin = (unsigned int)(tileY * blockSize);
							const unsigned int xEnd = std::min((unsigned int)((tileX + 1) * blockSize), width);
							const unsigned int yEnd = std::min((unsigned int)((tileY + 1) * blockSize), height);

							const int bi = tileId;

							for (unsigned int y = yBegin; y < yEnd; ++y)
							{
								for (unsigned int x = xBegin; x < xEnd; ++x)
								{
									Colour pixelAccum(0.0f, 0.0f, 0.0f);

									for (int s = 0; s < pilotSPP; ++s)
									{
										const float px = (float)x + sampler->next();
										const float py = (float)y + sampler->next();

										Ray ray = scene->camera.generateRay(px, py, sampler);
										Colour throughput(1.0f, 1.0f, 1.0f);
										Colour sampleCol = pathTrace(ray, throughput, 0, sampler);

										sampleCol.r = std::min(sampleCol.r, 100.0f);
										sampleCol.g = std::min(sampleCol.g, 100.0f);
										sampleCol.b = std::min(sampleCol.b, 100.0f);

										pixelAccum = pixelAccum + sampleCol;
									}

									Colour pilotMean = pixelAccum / (float)pilotSPP;
									const unsigned int idx = y * width + x;
									frameBuffer[idx] = pilotMean;

									const float lum = pilotMean.Lum();
									blockSumLum[bi] += lum;
									blockSumLumSq[bi] += lum * lum;
									blockPixelCount[bi]++;
								}
							}
						}
					});
			}

			for (size_t i = 0; i < workers.size(); ++i)
			{
				workers[i].join();
			}
		}

		// 2) Compute variance per block
		float maxVar = 0.0f;
		for (int bi = 0; bi < blockCount; bi++)
		{
			if (blockPixelCount[bi] <= 0)
			{
				blockVariance[bi] = 0.0f;
				continue;
			}

			const float invN = 1.0f / (float)blockPixelCount[bi];
			const float mean = blockSumLum[bi] * invN;
			const float meanSq = blockSumLumSq[bi] * invN;
			float var = meanSq - (mean * mean);
			if (var < 0.0f)
			{
				var = 0.0f;
			}

			blockVariance[bi] = var;
			if (var > maxVar)
			{
				maxVar = var;
			}
		}

		// 3) Allocate more samples to noisier blocks
		for (int bi = 0; bi < blockCount; bi++)
		{
			if (maxVar <= 1e-8f)
			{
				blockSPP[bi] = pilotSPP;
				continue;
			}

			float t = blockVariance[bi] / maxVar;
			int spp = pilotSPP + (int)floorf(t * (float)(maxSPPPerPixelThisFrame - pilotSPP) + 0.5f);
			if (spp < pilotSPP)
			{
				spp = pilotSPP;
			}
			if (spp > maxSPPPerPixelThisFrame)
			{
				spp = maxSPPPerPixelThisFrame;
			}
			blockSPP[bi] = spp;
		}
		//10-2
		// 4) Adaptive refinement pass
		{
			const int workerCount = std::max(1, numProcs);
			std::atomic<int> nextTile(0);
			std::vector<std::thread> workers;
			workers.reserve((size_t)workerCount);

			for (int threadIndex = 0; threadIndex < workerCount; ++threadIndex)
			{
				workers.emplace_back([&, threadIndex]()
					{
						Sampler* sampler = &samplers[threadIndex];

						while (true)
						{
							const int tileId = nextTile.fetch_add(1);
							if (tileId >= blockCount)
							{
								break;
							}

							const int tileX = tileId % blocksX;
							const int tileY = tileId / blocksX;

							const unsigned int xBegin = (unsigned int)(tileX * blockSize);
							const unsigned int yBegin = (unsigned int)(tileY * blockSize);
							const unsigned int xEnd = std::min((unsigned int)((tileX + 1) * blockSize), width);
							const unsigned int yEnd = std::min((unsigned int)((tileY + 1) * blockSize), height);

							const int bi = tileId;
							const int targetSPP = blockSPP[bi];
							if (targetSPP <= pilotSPP)
							{
								continue;
							}

							for (unsigned int y = yBegin; y < yEnd; ++y)
							{
								for (unsigned int x = xBegin; x < xEnd; ++x)
								{
									const unsigned int idx = y * width + x;
									Colour pixelAccum = frameBuffer[idx] * (float)pilotSPP;

									for (int s = pilotSPP; s < targetSPP; ++s)
									{
										const float px = (float)x + sampler->next();
										const float py = (float)y + sampler->next();

										Ray ray = scene->camera.generateRay(px, py, sampler);
										Colour throughput(1.0f, 1.0f, 1.0f);
										Colour sampleCol = pathTrace(ray, throughput, 0, sampler);

										sampleCol.r = std::min(sampleCol.r, 100.0f);
										sampleCol.g = std::min(sampleCol.g, 100.0f);
										sampleCol.b = std::min(sampleCol.b, 100.0f);

										pixelAccum = pixelAccum + sampleCol;
									}

									frameBuffer[idx] = pixelAccum / (float)targetSPP;
								}
							}
						}
					});
			}

			for (size_t i = 0; i < workers.size(); ++i)
			{
				workers[i].join();
			}
		}

		// 5) Commit to film + fill AOV/OIDN input buffers
		for (unsigned int y = 0; y < height; y++)
		{
			for (unsigned int x = 0; x < width; x++)
			{
				const float px = (float)x + 0.5f;
				const float py = (float)y + 0.5f;
				const unsigned int idx = y * width + x;
				const unsigned int aovBase = idx * 3;

				filmSplatMutex.lock();
				film->splat(px, py, frameBuffer[idx]);
				filmSplatMutex.unlock();

				Ray aovRay = scene->camera.generateRay(px, py, &samplers[0]);

				Colour alb = albedo(aovRay);
				albedoBuffer[aovBase + 0] = alb.r;
				albedoBuffer[aovBase + 1] = alb.g;
				albedoBuffer[aovBase + 2] = alb.b;

				Colour nor = viewNormals(aovRay);
				normalBuffer[aovBase + 0] = nor.r;
				normalBuffer[aovBase + 1] = nor.g;
				normalBuffer[aovBase + 2] = nor.b;

				const float invSPP = (film->SPP > 0) ? (1.0f / (float)film->SPP) : 0.0f;
				colorBuffer[aovBase + 0] = film->film[idx].r * invSPP;
				colorBuffer[aovBase + 1] = film->film[idx].g * invSPP;
				colorBuffer[aovBase + 2] = film->film[idx].b * invSPP;
			}
		}

		// 6) OIDN denoise: color + albedo + normal -> output
		{
			oidn::DeviceRef device = oidn::newDevice();
			device.commit();

			oidn::FilterRef filter = device.newFilter("RT");
			filter.setImage("color", colorBuffer.data(), oidn::Format::Float3, width, height);
			filter.setImage("albedo", albedoBuffer.data(), oidn::Format::Float3, width, height);
			filter.setImage("normal", normalBuffer.data(), oidn::Format::Float3, width, height);
			filter.setImage("output", outputBuffer.data(), oidn::Format::Float3, width, height);
			filter.set("hdr", true);
			filter.commit();
			filter.execute();

			const char* errorMessage = NULL;
			if (device.getError(errorMessage) != oidn::Error::None)
			{
				outputBuffer = colorBuffer;
			}
		}

		// 7) Tonemap denoised output and draw to screen
		for (unsigned int y = 0; y < height; y++)
		{
			for (unsigned int x = 0; x < width; x++)
			{
				const unsigned int idx = y * width + x;
				const unsigned int base = idx * 3;

				unsigned char r;
				unsigned char g;
				unsigned char b;
				tonemapLinearHDR(
					outputBuffer[base + 0],
					outputBuffer[base + 1],
					outputBuffer[base + 2],
					r, g, b, 1.0f);

				canvas->draw(x, y, r, g, b);
			}
		}
	}
	int getSPP()
	{
		return film->SPP;
	}
	void saveHDR(std::string filename)
	{
		film->save(filename);
	}
	void savePNG(std::string filename)
	{
		stbi_write_png(filename.c_str(), canvas->getWidth(), canvas->getHeight(), 3, canvas->getBackBuffer(), canvas->getWidth() * 3);
	}
	//8-2 给定场景中一个点 p、该点的法线 n、以及到达该点的颜色 col，判断该点是否能被相机看到，若能则将其贡献累积到 film 上。
	void connectToCamera(const Vec3& p, const Vec3& n, const Colour& col)
	{
		float x = 0.0f;
		float y = 0.0f;

		// 1) Project to camera film
		if (!scene->camera.projectOntoCamera(p, x, y))
		{
			return;
		}

		// 2) Visibility test
		if (!scene->visible(p, scene->camera.origin))
		{
			return;
		}

		// 3) Direction from point to camera
		Vec3 toCamera = scene->camera.origin - p;
		float distSq = toCamera.lengthSq();
		if (distSq <= (EPSILON * EPSILON))
		{
			return;
		}
		float dist = sqrtf(distSq);
		Vec3 wi = toCamera / dist;

		// 4) Geometry term + sensor importance
		float cosAtCamera = Dot(scene->camera.viewDirection, -wi);
		if (cosAtCamera <= 0.0f)
		{
			return;
		}

		float cosAtSurface = Dot(n, wi);
		if (cosAtSurface <= 0.0f)
		{
			return;
		}

		if (scene->camera.Afilm <= 0.0f)
		{
			return;
		}

		float cosAtCamera2 = cosAtCamera * cosAtCamera;
		float cosAtCamera4 = cosAtCamera2 * cosAtCamera2;

		float G = (cosAtCamera * cosAtSurface) / distSq;
		float We = 1.0f / (scene->camera.Afilm * cosAtCamera4);
		float weight = We * G;

		// 5) Accumulate to film
		Colour contribution = col * weight;
		film->splat(x, y, contribution);
	}
	//8-3 从光源出发，采样一个起始位置和方向，将起始点连接到相机，然后启动光路追踪。
	void lightTrace(Sampler* sampler, float scalePerPath)
	{
		// 1) Sample a light (same strategy as computeDirect)
		float lightPMF = 0.0f;
		Light* light = scene->sampleLight(sampler, lightPMF);
		if (light == NULL || lightPMF <= 0.0f)
		{
			return;
		}

		// 2) Sample position and direction from the light
		float pdfPosition = 0.0f;
		float pdfDirection = 0.0f;
		Vec3 p = light->samplePositionFromLight(sampler, pdfPosition);
		Vec3 wi = light->sampleDirectionFromLight(sampler, pdfDirection);

		if (pdfPosition <= 0.0f || pdfDirection <= 0.0f)
		{
			return;
		}

		// 3) Compute light-side normal
		ShadingData lightShadingData = {};
		lightShadingData.x = p;
		lightShadingData.wo = -wi;
		Vec3 lightNormal = light->normal(lightShadingData, wi);

		// 4) Initial radiance/colour weight
		Colour Le = light->evaluate(-wi);
		Colour col = Le / (pdfPosition * lightPMF);

		// 5) Connect starting point to camera (scaled per path)
		connectToCamera(p, lightNormal, col * scalePerPath);

		// 6) Build initial ray from light
		Ray ray;
		ray.init(p + (wi * EPSILON), wi);

		// 7) Initial throughput
		Colour pathThroughput = col / pdfDirection;

		// 8) Continue light path tracing
		lightTracePath(ray, pathThroughput, Le, sampler, scalePerPath);
	}
	// 8-4 递归追踪从光源出发的光路。与 pathTrace() 的主要区别是：不计算直接光照，而是在每个交点都尝试连接到相机。
	void lightTracePath(Ray& ray, Colour pathThroughput, Colour Le, Sampler* sampler, float scalePerPath, int depth = 0)
	{
		// Depth guard
		const int maxDepth = 64;
		if (depth > maxDepth)
		{
			return;
		}

		// 1) Intersect scene
		IntersectionData intersection = scene->traverse(ray);
		ShadingData shadingData = scene->calculateShadingData(intersection, ray);

		// 2) Miss => terminate
		if (shadingData.t >= FLT_MAX)
		{
			return;
		}

		// 3) Connect current hit point to camera
		Vec3 wi = scene->camera.origin - shadingData.x;
		float wiLenSq = wi.lengthSq();
		if (wiLenSq > (EPSILON * EPSILON))
		{
			wi = wi.normalize();
			Colour bsdfVal = shadingData.bsdf->evaluate(shadingData, wi);
			Colour col = pathThroughput * bsdfVal;
			connectToCamera(shadingData.x, shadingData.sNormal, col * scalePerPath);
		}

		// 4) Russian roulette
		float continueProb = std::min(0.95f, pathThroughput.Lum());
		if (continueProb <= 0.0f)
		{
			return;
		}
		if (sampler->next() > continueProb)
		{
			return;
		}
		pathThroughput = pathThroughput / continueProb;

		// 5) Sample next direction from BSDF
		Colour f;
		float pdf = 0.0f;
		Vec3 nextDir = shadingData.bsdf->sample(shadingData, sampler, f, pdf);
		if (pdf <= 0.0f)
		{
			return;
		}

		// 6) Update throughput
		float cosTheta = fabsf(Dot(shadingData.sNormal, nextDir));
		pathThroughput = pathThroughput * f * (cosTheta / pdf);

		// 7) Spawn next ray and recurse
		Ray nextRay;
		nextRay.init(shadingData.x + (nextDir * EPSILON), nextDir);
		lightTracePath(nextRay, pathThroughput, Le, sampler, scalePerPath, depth + 1);
	}
	// 9-1
	const std::vector<float>& getAlbedoBuffer() const
	{
		return albedoBuffer;
	}
	const std::vector<float>& getNormalBuffer() const
	{
		return normalBuffer;
	}

	// 9-2
	const std::vector<float>& getColorBuffer() const
	{
		return colorBuffer;
	}
	//9-3
	void tonemapLinearHDR(float inR, float inG, float inB, unsigned char& r, unsigned char& g, unsigned char& b, float exposure = 1.0f) const
	{
		inR = std::max(0.0f, inR * exposure);
		inG = std::max(0.0f, inG * exposure);
		inB = std::max(0.0f, inB * exposure);

		const float whitePoint = 4.0f;
		const float Cw = whitePoint / (1.0f + whitePoint);

		float outR = (inR / (1.0f + inR)) / Cw;
		float outG = (inG / (1.0f + inG)) / Cw;
		float outB = (inB / (1.0f + inB)) / Cw;

		const float invGamma = 1.0f / 2.2f;
		outR = powf(std::max(0.0f, outR), invGamma);
		outG = powf(std::max(0.0f, outG), invGamma);
		outB = powf(std::max(0.0f, outB), invGamma);

		outR = std::min(1.0f, std::max(0.0f, outR));
		outG = std::min(1.0f, std::max(0.0f, outG));
		outB = std::min(1.0f, std::max(0.0f, outB));

		r = (unsigned char)(outR * 255.0f);
		g = (unsigned char)(outG * 255.0f);
		b = (unsigned char)(outB * 255.0f);
	}
};

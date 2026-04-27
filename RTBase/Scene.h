#pragma once

#include "Core.h"
#include "Sampling.h"
#include "Geometry.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"
// 6-2 之前的代码里，虚拟相机是一个完美的“小孔成像（Pinhole）”相机，不管远近，所有的东西拍出来都极其清晰。但这不真实。真实世界的高级单反相机是有**镜头（Lens）的，能拍出景深（Depth of Field）**效果——也就是“对焦的地方清晰，背景或前景模糊”。这一步就是要模拟这个物理镜头
class Camera
{
public:
	Matrix projectionMatrix;
	Matrix inverseProjectionMatrix;
	Matrix camera;
	Matrix cameraToView;
	float width = 0;
	float height = 0;
	Vec3 origin;
	Vec3 viewDirection;
	float Afilm;

	// Thin-lens parameters
	float lensRadius = 0.0f;   // 0 => pinhole fallback
	float focalDistance = 5.0f;

	void init(Matrix ProjectionMatrix, int screenwidth, int screenheight)
	{
		projectionMatrix = ProjectionMatrix;
		inverseProjectionMatrix = ProjectionMatrix.invert();
		width = (float)screenwidth;
		height = (float)screenheight;
		float Wlens = (2.0f / ProjectionMatrix.a[1][1]);
		float aspect = ProjectionMatrix.a[0][0] / ProjectionMatrix.a[1][1];
		float Hlens = Wlens * aspect;
		Afilm = Wlens * Hlens;
	}
	void updateView(Matrix V)
	{
		camera = V;
		cameraToView = V.invert();
		origin = camera.mulPoint(Vec3(0, 0, 0));
		viewDirection = inverseProjectionMatrix.mulPointAndPerspectiveDivide(Vec3(0, 0, 1));
		viewDirection = camera.mulVec(viewDirection);
		viewDirection = viewDirection.normalize();
	}

	void setThinLens(float _lensRadius, float _focalDistance)
	{
		lensRadius = std::max(0.0f, _lensRadius);
		focalDistance = std::max(EPSILON, _focalDistance);
	}
	// Add code here:1-1 这是光线追踪的第一步。你需要把屏幕上的一个二维像素点（Pixel），反向转换成3D虚拟世界里的一条真实光线（Ray）
	Ray generateRay(float x, float y, Sampler* sampler)
	{
		// Screen space -> NDC [-1, 1]
		float ndcX = (2.0f * x / width) - 1.0f;
		float ndcY = 1.0f - (2.0f * y / height);

		// NDC -> view space
		Vec3 viewPos = inverseProjectionMatrix.mulPointAndPerspectiveDivide(Vec3(ndcX, ndcY, 1.0f));

		// Center ray direction (pinhole direction)
		Vec3 pinholeDir = camera.mulVec(viewPos).normalize();

		// Fallback to pinhole when DOF disabled
		if (sampler == NULL || lensRadius <= 0.0f)
		{
			return Ray(origin, pinholeDir);
		}

		// Compute focus point on focal plane:
		// plane normal = viewDirection, distance = focalDistance from camera origin
		float denom = Dot(pinholeDir, viewDirection);
		if (fabsf(denom) < 1e-6f)
		{
			return Ray(origin, pinholeDir);
		}
		float tFocus = focalDistance / denom;
		if (tFocus <= EPSILON)
		{
			tFocus = focalDistance;
		}
		Vec3 xFocus = origin + (pinholeDir * tFocus);

		// Sample point on lens disk (uniform area)
		float u1 = sampler->next();
		float u2 = sampler->next();
		float r = lensRadius * sqrtf(u1);
		float phi = 2.0f * (float)M_PI * u2;

		Vec3 lensLocal(r * cosf(phi), r * sinf(phi), 0.0f);
		Vec3 lensWorld = origin + camera.mulVec(lensLocal);

		Vec3 dir = (xFocus - lensWorld).normalize();
		return Ray(lensWorld, dir);
	}

	// Keep old call-site compatibility
	Ray generateRay(float x, float y)
	{
		return generateRay(x, y, NULL);
	}

	bool projectOntoCamera(const Vec3& p, float& x, float& y)
	{
		Vec3 pview = cameraToView.mulPoint(p);
		Vec3 pproj = projectionMatrix.mulPointAndPerspectiveDivide(pview);
		x = (pproj.x + 1.0f) * 0.5f;
		y = (pproj.y + 1.0f) * 0.5f;
		if (x < 0 || x > 1.0f || y < 0 || y > 1.0f)
		{
			return false;
		}
		x = x * width;
		y = 1.0f - y;
		y = y * height;
		return true;
	}
};

class Scene
{
public:
	std::vector<Triangle> triangles;
	std::vector<BSDF*> materials;
	std::vector<Light*> lights;
	Light* background = NULL;
	BVHNode* bvh = NULL;
	Camera camera;
	AABB bounds;
	void build()
	{
		// 3-2有了一堆杂乱无章的三角形和空架子，现在要把三角形装进盒子里，并且不断地把大盒子切分成小盒子。这里最关键的是**“怎么切”才最快**？讲义要求你实现一种叫 SAH（表面积启发式算法，Surface Area Heuristic） 的高级策略。它的核心思想是：盒子的表面积越大，越容易被光线射中，所以我们要根据面积和三角形的数量来精打细算，找到最完美的切分位置
		// Build BVH once at scene startup
		if (!triangles.empty())
		{
			bvh = new BVHNode();
			bvh->build(triangles);
		}

		// Do not touch the code below this line!
		// Build light list
		for (int i = 0; i < triangles.size(); i++)
		{
			if (materials[triangles[i].materialIndex]->isLight())
			{
				AreaLight* light = new AreaLight();
				light->triangle = &triangles[i];
				light->emission = materials[triangles[i].materialIndex]->emission;
				lights.push_back(light);
			}
		}
	}
	IntersectionData traverse(const Ray& ray)
	{
		// 3-3 
		// 优先走 BVH
		if (bvh != NULL)
		{
			return bvh->traverse(ray, triangles);
		}

		// fallback：线性遍历
		IntersectionData intersection;
		intersection.t = FLT_MAX;
		for (int i = 0; i < triangles.size(); i++)
		{
			float t;
			float u;
			float v;
			if (triangles[i].rayIntersect(ray, t, u, v))
			{
				if (t < intersection.t)
				{
					intersection.t = t;
					intersection.ID = i;
					intersection.alpha = u;
					intersection.beta = v;
					intersection.gamma = 1.0f - (u + v);
				}
			}
		}
		return intersection;
	}
	// 5-1  如果场景里有几十上百盏灯，为了算出物体受到的光照，如果每次光线碰撞都要把所有的灯算一遍，电脑会卡死。为了提速，我们需要用到前面讲过的蒙特卡洛“扔骰子”的方法——每次只随机挑选一盏灯来探测
	Light* sampleLight(Sampler* sampler, float& pmf)
	{
		const int lightCount = (int)lights.size();
		if (lightCount == 0)
		{
			pmf = 0.0f;
			return NULL;
		}

		const float xi = sampler->next();
		int lightIndex = (int)floorf(xi * (float)lightCount);

		// Safety clamp (in case xi hits 1.0 due to numeric edge cases)
		if (lightIndex >= lightCount)
		{
			lightIndex = lightCount - 1;
		}
		if (lightIndex < 0)
		{
			lightIndex = 0;
		}

		pmf = 1.0f / (float)lightCount;
		return lights[lightIndex];
	}
	// Do not modify any code below this line
	void init(std::vector<Triangle> meshTriangles, std::vector<BSDF*> meshMaterials, Light* _background)
	{
		for (int i = 0; i < meshTriangles.size(); i++)
		{
			triangles.push_back(meshTriangles[i]);
			bounds.extend(meshTriangles[i].vertices[0].p);
			bounds.extend(meshTriangles[i].vertices[1].p);
			bounds.extend(meshTriangles[i].vertices[2].p);
		}
		for (int i = 0; i < meshMaterials.size(); i++)
		{
			materials.push_back(meshMaterials[i]);
		}
		background = _background;
		if (background->totalIntegratedPower() > 0)
		{
			lights.push_back(background);
		}
	}
	// 3-4
	bool visible(const Vec3& p1, const Vec3& p2)
	{
		Ray ray;
		Vec3 dir = p2 - p1;
		float maxT = dir.length() - (2.0f * EPSILON);
		if (maxT <= EPSILON)
		{
			return true;
		}
		dir = dir.normalize();
		ray.init(p1 + (dir * EPSILON), dir);

		// Use BVH when available
		if (bvh != NULL)
		{
			return bvh->traverseVisible(ray, triangles, maxT);
		}

		// Fallback: linear visibility test
		for (int i = 0; i < triangles.size(); i++)
		{
			float t;
			float u;
			float v;
			if (triangles[i].rayIntersect(ray, t, u, v))
			{
				if (t > EPSILON && t < maxT)
				{
					return false;
				}
			}
		}
		return true;
	}
	Colour emit(Triangle* light, ShadingData shadingData, Vec3 wi)
	{
		return materials[light->materialIndex]->emit(shadingData, wi);
	}
	ShadingData calculateShadingData(IntersectionData intersection, Ray& ray)
	{
		ShadingData shadingData = {};
		if (intersection.t < FLT_MAX)
		{
			shadingData.x = ray.at(intersection.t);
			shadingData.gNormal = triangles[intersection.ID].gNormal();
			triangles[intersection.ID].interpolateAttributes(intersection.alpha, intersection.beta, intersection.gamma, shadingData.sNormal, shadingData.tu, shadingData.tv);
			shadingData.bsdf = materials[triangles[intersection.ID].materialIndex];
			shadingData.wo = -ray.dir;
			if (shadingData.bsdf->isTwoSided())
			{
				if (Dot(shadingData.wo, shadingData.sNormal) < 0)
				{
					shadingData.sNormal = -shadingData.sNormal;
				}
				if (Dot(shadingData.wo, shadingData.gNormal) < 0)
				{
					shadingData.gNormal = -shadingData.gNormal;
				}
			}
			shadingData.frame.fromVector(shadingData.sNormal);
			shadingData.t = intersection.t;
		} else
		{
			shadingData.wo = -ray.dir;
			shadingData.t = intersection.t;
		}
		return shadingData;
	}
};
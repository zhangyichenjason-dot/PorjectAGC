#pragma once

#include "Core.h"
#include "Geometry.h"
#include "Materials.h"
#include "Sampling.h"

#pragma warning( disable : 4244)

class SceneBounds
{
public:
	Vec3 sceneCentre;
	float sceneRadius;
};

class Light
{
public:
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf) = 0;
	virtual Colour evaluate(const Vec3& wi) = 0;
	virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual bool isArea() = 0;
	virtual Vec3 normal(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float totalIntegratedPower() = 0;
	virtual Vec3 samplePositionFromLight(Sampler* sampler, float& pdf) = 0;
	virtual Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf) = 0;
};

class AreaLight : public Light
{
public:
	Triangle* triangle = NULL;
	Colour emission;
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf)
	{
		emittedColour = emission;
		return triangle->sample(sampler, pdf);
	}
	Colour evaluate(const Vec3& wi)
	{
		if (Dot(wi, triangle->gNormal()) < 0)
		{
			return emission;
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		return 1.0f / triangle->area;
	}
	bool isArea()
	{
		return true;
	}
	Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
	{
		return triangle->gNormal();
	}
	float totalIntegratedPower()
	{
		return (triangle->area * emission.Lum());
	}
	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		return triangle->sample(sampler, pdf);
	}
	// 8-1（RenderingAlgorithms）要实现光追。从面光源法线方向的半球上进行余弦加权采样（cosine sampling）。采样方向应在面光源几何法源所定义的局部坐标系中生成，然后转换回世界坐标。
	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	{
		Frame frame;
		frame.fromVector(triangle->gNormal());

		float r1 = sampler->next();
		float r2 = sampler->next();

		Vec3 wiLocal = SamplingDistributions::cosineSampleHemisphere(r1, r2);
		pdf = SamplingDistributions::cosineHemispherePDF(wiLocal);

		return frame.toWorld(wiLocal);
	}
};

class BackgroundColour : public Light
{
public:
	Colour emission;
	BackgroundColour(Colour _emission)
	{
		emission = _emission;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		pdf = SamplingDistributions::uniformSpherePDF(wi);
		reflectedColour = emission;
		return wi;
	}
	Colour evaluate(const Vec3& wi)
	{
		return emission;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		return SamplingDistributions::uniformSpherePDF(wi);
	}
	bool isArea()
	{
		return false;
	}
	Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
	{
		return -wi;
	}
	float totalIntegratedPower()
	{
		return emission.Lum() * 4.0f * M_PI;
	}
	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		Vec3 p = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		p = p * use<SceneBounds>().sceneRadius;
		p = p + use<SceneBounds>().sceneCentre;
		pdf = 4 * M_PI * use<SceneBounds>().sceneRadius * use<SceneBounds>().sceneRadius;
		return p;
	}
	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	{
		Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		pdf = SamplingDistributions::uniformSpherePDF(wi);
		return wi;
	}
};

class EnvironmentMap : public Light
{
public:
	Texture* env;

	//11-1(作业要求：环境贴图重要性采样 + MIS )存储预计算好的 CDF 表，供采样时查找用
	std::vector<float> marginalCDF;                    // size: H + 1
	std::vector<float> marginalPDF;                    // size: H
	std::vector<std::vector<float>> conditionalCDF;    // H rows, each size: W + 1
	std::vector<std::vector<float>> conditionalPDF;    // H rows, each size: W
	float totalWeight = 0.0f;                          // sum of all pixel weights

	//11-2 把环境贴图的亮度分布预先转成可以快速查找的累积概率表，只算一次
	EnvironmentMap(Texture* _env)
	{
		env = _env;

		const int W = env->width;
		const int H = env->height;

		marginalCDF.assign(H + 1, 0.0f);
		marginalPDF.assign(H, 0.0f);
		conditionalCDF.assign(H, std::vector<float>(W + 1, 0.0f));
		conditionalPDF.assign(H, std::vector<float>(W, 0.0f));
		totalWeight = 0.0f;

		// Step 1 + Step 2: per-row weights, conditional CDF/PDF over columns
		for (int j = 0; j < H; j++)
		{
			const float sinTheta = sinf((((float)j + 0.5f) / (float)H) * M_PI);

			float rowSum = 0.0f;
			for (int i = 0; i < W; i++)
			{
				const float w = env->texels[(j * W) + i].Lum() * sinTheta;
				rowSum += w;
			}

			conditionalCDF[j][0] = 0.0f;

			if (rowSum > 0.0f)
			{
				float accum = 0.0f;
				for (int i = 0; i < W; i++)
				{
					const float w = env->texels[(j * W) + i].Lum() * sinTheta;
					const float pDiscrete = w / rowSum; // discrete prob for pixel i in row j
					accum += pDiscrete;

					conditionalCDF[j][i + 1] = accum;                  // [0,1] CDF
					conditionalPDF[j][i] = pDiscrete * (float)W;       // continuous PDF in u
				}
				conditionalCDF[j][W] = 1.0f;
			}
			else
			{
				// Fallback: uniform along row
				for (int i = 0; i < W; i++)
				{
					conditionalCDF[j][i + 1] = (float)(i + 1) / (float)W;
					conditionalPDF[j][i] = 1.0f;
				}
			}

			marginalPDF[j] = rowSum; // temporary: store row weight
			totalWeight += rowSum;
		}

		// Step 3: marginal CDF/PDF over rows
		marginalCDF[0] = 0.0f;
		if (totalWeight > 0.0f)
		{
			float accum = 0.0f;
			for (int j = 0; j < H; j++)
			{
				const float pDiscrete = marginalPDF[j] / totalWeight; // discrete prob for row j
				accum += pDiscrete;

				marginalCDF[j + 1] = accum;                // [0,1] CDF
				marginalPDF[j] = pDiscrete * (float)H;     // continuous PDF in v
			}
			marginalCDF[H] = 1.0f;
		}
		else
		{
			// Fallback: uniform rows
			for (int j = 0; j < H; j++)
			{
				marginalCDF[j + 1] = (float)(j + 1) / (float)H;
				marginalPDF[j] = 1.0f;
			}
			totalWeight = 1.0f;
		}
	}
	//11-3 用两个随机数，通过查 CDF 表找到一个按亮度加权的方向
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		const int W = env->width;
		const int H = env->height;

		// Safety fallback
		if (W <= 0 || H <= 0 || marginalCDF.empty() || conditionalCDF.empty())
		{
			Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
			pdf = SamplingDistributions::uniformSpherePDF(wi);
			reflectedColour = evaluate(wi);
			return wi;
		}

		const float r1 = sampler->next();
		const float r2 = sampler->next();

		// 1) Sample row j from marginal CDF
		int low = 0;
		int high = (int)marginalCDF.size() - 1; // H
		while (low + 1 < high)
		{
			const int mid = (low + high) / 2;
			if (marginalCDF[mid] <= r1)
			{
				low = mid;
			}
			else
			{
				high = mid;
			}
		}
		int j = low;
		if (j < 0)
		{
			j = 0;
		}
		if (j >= H)
		{
			j = H - 1;
		}

		// 2) Sample column i from conditional CDF of row j
		const std::vector<float>& rowCDF = conditionalCDF[j];
		low = 0;
		high = (int)rowCDF.size() - 1; // W
		while (low + 1 < high)
		{
			const int mid = (low + high) / 2;
			if (rowCDF[mid] <= r2)
			{
				low = mid;
			}
			else
			{
				high = mid;
			}
		}
		int i = low;
		if (i < 0)
		{
			i = 0;
		}
		if (i >= W)
		{
			i = W - 1;
		}

		// 3) Pixel center -> spherical direction
		const float phi = (((float)i + 0.5f) / (float)W) * 2.0f * M_PI;
		const float theta = (((float)j + 0.5f) / (float)H) * M_PI;

		const float sinTheta = sinf(theta);
		const float cosTheta = cosf(theta);
		const float cosPhi = cosf(phi);
		const float sinPhi = sinf(phi);

		Vec3 wi(cosPhi * sinTheta, cosTheta, sinPhi * sinTheta);

		// 4) Convert (u,v) continuous PDF to solid-angle PDF
		if (sinTheta <= 1e-6f)
		{
			pdf = 0.0f;
		}
		else
		{
			const float uvPdf = marginalPDF[j] * conditionalPDF[j][i];
			pdf = uvPdf / (2.0f * M_PI * M_PI * sinTheta);
		}

		reflectedColour = evaluate(wi);
		return wi;
	}
	Colour evaluate(const Vec3& wi)
	{
		float u = atan2f(wi.z, wi.x);
		u = (u < 0.0f) ? u + (2.0f * M_PI) : u;
		u = u / (2.0f * M_PI);
		float v = acosf(wi.y) / M_PI;
		return env->sample(u, v);
	}
	//11-4给定任意一个方向，返回它被重要性采样采到的概率密度，MIS 时需要用这个值
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		const int W = env->width;
		const int H = env->height;

		if (W <= 0 || H <= 0 || (int)marginalPDF.size() != H || (int)conditionalPDF.size() != H)
		{
			return SamplingDistributions::uniformSpherePDF(wi);
		}

		float y = wi.y;
		if (y < -1.0f)
		{
			y = -1.0f;
		}
		else if (y > 1.0f)
		{
			y = 1.0f;
		}

		float phi = atan2f(wi.z, wi.x);
		if (phi < 0.0f)
		{
			phi += 2.0f * M_PI;
		}
		const float u = phi / (2.0f * M_PI);

		const float theta = acosf(y);
		const float v = theta / M_PI;

		int i = (int)(u * (float)W);
		int j = (int)(v * (float)H);

		if (i < 0)
		{
			i = 0;
		}
		else if (i >= W)
		{
			i = W - 1;
		}

		if (j < 0)
		{
			j = 0;
		}
		else if (j >= H)
		{
			j = H - 1;
		}

		const float sinTheta = sinf(theta);
		if (sinTheta <= 1e-6f)
		{
			return 0.0f;
		}

		const float uvPdf = marginalPDF[j] * conditionalPDF[j][i];
		return uvPdf / (2.0f * M_PI * M_PI * sinTheta);
	}
	bool isArea()
	{
		return false;
	}
	Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
	{
		return -wi;
	}
	float totalIntegratedPower()
	{
		float total = 0;
		for (int i = 0; i < env->height; i++)
		{
			float st = sinf(((float)i / (float)env->height) * M_PI);
			for (int n = 0; n < env->width; n++)
			{
				total += (env->texels[(i * env->width) + n].Lum() * st);
			}
		}
		total = total / (float)(env->width * env->height);
		return total * 4.0f * M_PI;
	}
	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		// Samples a point on the bounding sphere of the scene. Feel free to improve this.
		Vec3 p = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		p = p * use<SceneBounds>().sceneRadius;
		p = p + use<SceneBounds>().sceneCentre;
		pdf = 1.0f / (4 * M_PI * SQ(use<SceneBounds>().sceneRadius));
		return p;
	}
	//12-1（作业要求：环境贴图重要性采样 + MIS ）
	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	{
		const int W = env->width;
		const int H = env->height;

		// Safety fallback
		if (W <= 0 || H <= 0 || marginalCDF.empty() || conditionalCDF.empty())
		{
			Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
			pdf = SamplingDistributions::uniformSpherePDF(wi);
			return -wi;
		}

		const float r1 = sampler->next();
		const float r2 = sampler->next();

		// 1) Sample row j from marginal CDF
		int low = 0;
		int high = (int)marginalCDF.size() - 1; // H
		while (low + 1 < high)
		{
			const int mid = (low + high) / 2;
			if (marginalCDF[mid] <= r1)
			{
				low = mid;
			}
			else
			{
				high = mid;
			}
		}
		int j = low;
		if (j < 0)
		{
			j = 0;
		}
		if (j >= H)
		{
			j = H - 1;
		}

		// 2) Sample column i from conditional CDF of row j
		const std::vector<float>& rowCDF = conditionalCDF[j];
		low = 0;
		high = (int)rowCDF.size() - 1; // W
		while (low + 1 < high)
		{
			const int mid = (low + high) / 2;
			if (rowCDF[mid] <= r2)
			{
				low = mid;
			}
			else
			{
				high = mid;
			}
		}
		int i = low;
		if (i < 0)
		{
			i = 0;
		}
		if (i >= W)
		{
			i = W - 1;
		}

		// 3) Pixel center -> spherical direction (same as sample())
		const float phi = (((float)i + 0.5f) / (float)W) * 2.0f * M_PI;
		const float theta = (((float)j + 0.5f) / (float)H) * M_PI;

		const float sinTheta = sinf(theta);
		const float cosTheta = cosf(theta);
		const float cosPhi = cosf(phi);
		const float sinPhi = sinf(phi);

		Vec3 wi(cosPhi * sinTheta, cosTheta, sinPhi * sinTheta);

		// 4) Convert (u,v) continuous PDF to solid-angle PDF
		if (sinTheta <= 1e-6f)
		{
			pdf = 0.0f;
		}
		else
		{
			const float uvPdf = marginalPDF[j] * conditionalPDF[j][i];
			pdf = uvPdf / (2.0f * M_PI * M_PI * sinTheta);
		}

		// 5) For light tracing from env -> shoot into scene
		return -wi;
	}
};
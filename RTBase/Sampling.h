#pragma once

#include "Core.h"
#include <random>
#include <algorithm>

class Sampler
{
public:
	virtual ~Sampler() = default;
	virtual float next() = 0;
};
// 4-4 蒙特卡洛渲染的命脉就是“随机”。电脑自带的普通随机数生成器（比如线性的 LCG）其实存在隐藏的规律，这种规律会导致渲染出的画面上出现奇怪的条纹或网格。为了解决这个问题，你需要实现一个图形学界公认的高质量伪随机数生成器：Mersenne Twister（梅森旋转算法）
class MTRandom : public Sampler
{
public:
	std::mt19937 generator;
	std::uniform_real_distribution<float> dist;

	explicit MTRandom(unsigned int seed = 1u)
		: generator(seed), dist(0.0f, 1.0f)
	{
	}

	float next() override
	{
		return dist(generator);
	}
};
// 4-1  当一束光打到粗糙的墙面上时，它会向四面八方漫无目的地散开。为了模拟这个现象，我们需要在代码里生成随机方向的光线，让它们像刺猬的刺一样，在一个**半球（Hemisphere，代表物体表面向外）或完整的球体（Sphere）上均匀地（Uniformly）**分布
// Note all of these distributions assume z-up coordinate system
class SamplingDistributions
{
public:
	static Vec3 uniformSampleHemisphere(float r1, float r2)
	{
		// Uniform on hemisphere (z >= 0):
		// z = r1, phi = 2*pi*r2
		float z = r1;
		float phi = 2.0f * (float)M_PI * r2;
		float sinTheta = sqrtf(std::max(0.0f, 1.0f - (z * z)));
		float x = cosf(phi) * sinTheta;
		float y = sinf(phi) * sinTheta;
		return Vec3(x, y, z);
	}
	static float uniformHemispherePDF(const Vec3 wi)
	{
		// Constant over hemisphere
		return (wi.z >= 0.0f) ? (1.0f / (2.0f * (float)M_PI)) : 0.0f;
	}
	// 4-2  这是上一步的“进阶升级版”。在物理世界中，垂直打在表面上的光比倾斜擦过表面的光要强得多（余弦定理）。所以，我们在“扔骰子”生成随机光线时，不应该傻傻地均匀发射，而应该把更多的光线集中在垂直向上的方向，边缘方向少发一点，这就叫余弦采样
	static Vec3 cosineSampleHemisphere(float r1, float r2)
	{
		// Malley's method:
		// 1) Uniform sample on unit disk
		// 2) Project up to hemisphere
		float r = sqrtf(r1);
		float phi = 2.0f * (float)M_PI * r2;

		float x = r * cosf(phi);
		float y = r * sinf(phi);
		float z = sqrtf(std::max(0.0f, 1.0f - (x * x) - (y * y)));

		return Vec3(x, y, z);
	}
	static float cosineHemispherePDF(const Vec3 wi)
	{
		// Add code here
		return 1.0f;
	}
	static Vec3 uniformSampleSphere(float r1, float r2)
	{
		// Uniform on full sphere:
		// z = 1 - 2*r1, phi = 2*pi*r2
		float z = 1.0f - (2.0f * r1);
		float phi = 2.0f * (float)M_PI * r2;
		float sinTheta = sqrtf(std::max(0.0f, 1.0f - (z * z)));
		float x = cosf(phi) * sinTheta;
		float y = sinf(phi) * sinTheta;
		return Vec3(x, y, z);
	}
	static float uniformSpherePDF(const Vec3& wi)
	{
		(void)wi;
		// Constant over sphere
		return 1.0f / (4.0f * (float)M_PI);
	}
};
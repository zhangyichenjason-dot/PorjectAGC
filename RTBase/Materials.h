#pragma once

#include "Core.h"
#include "Imaging.h"
#include "Sampling.h"

#pragma warning( disable : 4244)
#pragma warning( disable : 4305) // Double to float

class BSDF;

class ShadingData
{
public:
	Vec3 x;
	Vec3 wo;
	Vec3 sNormal;
	Vec3 gNormal;
	float tu;
	float tv;
	Frame frame;
	BSDF* bsdf;
	float t;
	ShadingData() {}
	ShadingData(Vec3 _x, Vec3 n)
	{
		x = _x;
		gNormal = n;
		sNormal = n;
		bsdf = NULL;
	}
};

class ShadingHelper
{
public:
	// 7-4 金属（导体）的折射率是复数，虚部k描述光在金属内部的衰减（趋肤效应）。Rs（垂直偏振）和Rp（平行偏振）：光是横波，有两个偏振方向，金属对它们的反射率不同。实际渲染用非偏振光，取平均 对rgb三通道分别算：金属的颜色正是因为不同波长的光反射率不同（比如金反射红色多，蓝色少），所以必须分通道计算
	static Colour fresnelConductor(float cosTheta, Colour ior, Colour k)
	{
		cosTheta = std::max(0.0f, std::min(1.0f, cosTheta));

		float cosThetaSq = cosTheta * cosTheta;
		float sinThetaSq = 1.0f - cosThetaSq;
		(void)sinThetaSq;

		auto evalChannel = [&](float n, float kk)
			{
				float n2k2 = n * n + kk * kk;

				float Rs_num = n2k2 - 2.0f * n * cosTheta + cosThetaSq;
				float Rs_den = n2k2 + 2.0f * n * cosTheta + cosThetaSq;
				float Rs = Rs_num / Rs_den;

				float Rp_num = n2k2 * cosThetaSq - 2.0f * n * cosTheta + 1.0f;
				float Rp_den = n2k2 * cosThetaSq + 2.0f * n * cosTheta + 1.0f;
				float Rp = Rp_num / Rp_den;

				return (Rs + Rp) * 0.5f;
			};

		return Colour(
			evalChannel(ior.r, k.r),
			evalChannel(ior.g, k.g),
			evalChannel(ior.b, k.b)
		);
	}
	// 7-5 用于玻璃材质
	static float fresnelDielectric(float cosTheta, float iorInt, float iorExt)
	{
		cosTheta = std::max(-1.0f, std::min(1.0f, cosTheta));

		float etaI = iorExt;
		float etaT = iorInt;

		if (cosTheta < 0.0f)
		{
			cosTheta = -cosTheta;
			std::swap(etaI, etaT);
		}

		float sinThetaI = sqrtf(std::max(0.0f, 1.0f - cosTheta * cosTheta));
		float sinThetaT = (etaI / etaT) * sinThetaI;

		if (sinThetaT >= 1.0f)
		{
			return 1.0f;
		}

		float cosThetaT = sqrtf(std::max(0.0f, 1.0f - sinThetaT * sinThetaT));

		float Rs = ((etaI * cosTheta) - (etaT * cosThetaT)) /
			((etaI * cosTheta) + (etaT * cosThetaT));
		float Rp = ((etaT * cosTheta) - (etaI * cosThetaT)) /
			((etaT * cosTheta) + (etaI * cosThetaT));

		return 0.5f * (Rs * Rs + Rp * Rp);
	}
	// 7-6 Lambda函数是Smith遮蔽函数的核心，描述"在这个方向上，有多少比例的微表面被其他微表面挡住了"。先算tanThetaSq：因为GGX的遮蔽函数是关于tan(θ)的函数 sqrt(1 + alpha ^ 2 * tan2)：alpha是粗糙度，粗糙度越大，遮蔽越严重
	static float lambdaGGX(Vec3 wi, float alpha)
	{
		float cosTheta = fabsf(wi.z);
		if (cosTheta >= 1.0f)
		{
			return 0.0f;
		}

		float tanThetaSq = (1.0f - (cosTheta * cosTheta)) / (cosTheta * cosTheta);
		return (-1.0f + sqrtf(1.0f + (alpha * alpha * tanThetaSq))) * 0.5f;
	}
	// 7-7 Smith近似假设入射方向的遮蔽和出射方向的遮蔽统计独立，所以总遮蔽=两个方向各自遮蔽相乘。对入射方向wi和出射方向wo分别算Lambda，然后组合。G越接近1说明遮蔽越少（表面越光滑）。
	static float Gggx(Vec3 wi, Vec3 wo, float alpha)
	{
		return 1.0f / ((1.0f + lambdaGGX(wi, alpha)) * (1.0f + lambdaGGX(wo, alpha)));
	}
	// 7-8 NDF（法线分布函数）描述微表面中有多少比例的微法线朝向h方向。
	// h.z是微法线与宏观法线的夹角的余弦，越接近1说明h越接近法线，光滑表面几乎所有微法线都朝这个方向
	// alphaSq在分子：粗糙度越大，NDF越宽（能量分散到更多方向）
	// cos4(θ)在分母：是从投影面积到立体角的Jacobian
	static float Dggx(Vec3 h, float alpha)
	{
		float cosTheta = h.z;
		if (cosTheta <= 0.0f)
		{
			return 0.0f;
		}

		float cos2Theta = cosTheta * cosTheta;
		float tan2Theta = (1.0f - cos2Theta) / cos2Theta;
		float alphaSq = alpha * alpha;
		float denom = (float)M_PI * cos2Theta * cos2Theta * (alphaSq + tan2Theta) * (alphaSq + tan2Theta);

		return alphaSq / denom;
	}
};

class BSDF
{
public:
	Colour emission;
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) = 0;
	virtual Colour evaluate(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual bool isPureSpecular() = 0;
	virtual bool isTwoSided() = 0;
	bool isLight()
	{
		return emission.Lum() > 0 ? true : false;
	}
	void addLight(Colour _emission)
	{
		emission = _emission;
	}
	Colour emit(const ShadingData& shadingData, const Vec3& wi)
	{
		return emission;
	}
	virtual float mask(const ShadingData& shadingData) = 0;
};


class DiffuseBSDF : public BSDF
{
public:
	Texture* albedo;
	DiffuseBSDF() = default;
	DiffuseBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	// 7-2 余弦采样：不是均匀地往所有方向发射光线，而是偏向法线方向多发，边缘少发。这样采样分布和 cosθ 项匹配，减少方差。转回世界坐标：采样是在局部坐标（z朝上）做的，渲染器里用的是世界坐标，所以要转换
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = SamplingDistributions::cosineHemispherePDF(wi);
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	// 7-2 转到局部坐标：PDF公式在局部坐标下定义，传入的wi是世界坐标，要先转换才能正确计算
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class MirrorBSDF : public BSDF
{
public:
	Texture* albedo;
	MirrorBSDF() = default;
	MirrorBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	// 7-3 x、y取反，z不变：这是局部坐标下镜面反射的几何意义。法线朝z轴，反射就是把水平分量翻转，垂直分量（沿法线）保持不变。除以cosθ（即wiLocal.z）：Mirror的BSDF在数学上包含一个狄拉克delta函数，当与渲染方程的cosθ项和pdf组合后，最终留下来需要除掉的正是这个cosθ。这是mirror材质的标准处理
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		(void)sampler;
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		pdf = 1.0f;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) ;
		return shadingData.frame.toWorld(wiLocal);
	}
	// 7-3 镜面只在一个精确方向有响应，evaluate用于任意方向，永远为0
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		(void)shadingData;
		(void)wi;
		return Colour(0.0f, 0.0f, 0.0f);
	}
	// 7-3  同理，除了那个唯一方向外，概率密度为0
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		(void)shadingData;
		(void)wi;
		return 0.0f;
	}
	bool isPureSpecular()
	{
		return true;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class ConductorBSDF : public BSDF
{
public:
	Texture* albedo;
	Colour eta;
	Colour k;
	float alpha;
	ConductorBSDF() = default;
	ConductorBSDF(Texture* _albedo, Colour _eta, Colour _k, float roughness)
	{
		albedo = _albedo;
		eta = _eta;
		k = _k;
		alpha = 1.62142f * sqrtf(roughness);
	}
	// 7-9 evaluate（Cook-Torrance BRDF）：算半程向量h：微表面理论中，只有法线恰好朝向wi和wo中间（即h方向）的微表面才能把光从wo方向反射到wi方向
		// F：这个微表面的菲涅尔反射率，决定多少光被反射
		// D：有多少比例的微表面法线朝向h，D越大说明这个方向的微表面越多
		// G：考虑遮蔽，不是所有朝向h的微表面都可见
		// 除以(4 cosθi cosθo)：是从微表面局部坐标转换到宏观渲染方程的Jacobian

		// sample：先采样微表面法线h，而不是直接采样wi：因为BRDF峰值在h集中的地方，按NDF采样h是重要性采样cosTheta_h的公式：是GGX分布的逆CDF，把均匀随机数r1映射到符合GGX分布的角度对h做反射得到wi：给定微表面法线h和出射方向wo，镜面反射唯一确定wipdf = Dh.z / (4dot(wo, h))：这是从h的分布转换到wi的分布的Jacobian，因为同样一个微法线采样，可能对应不同的立体角alpha < epsilon退化为Mirror： 完全光滑时GGX退化为狄拉克delta，直接走镜面处理更稳定
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);

		if (alpha < 0.001f)
		{
			Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			float cosTheta = fabsf(wiLocal.z);
			if (cosTheta <= 0.0f)
			{
				pdf = 0.0f;
				reflectedColour = Colour(0.0f, 0.0f, 0.0f);
				return Vec3(0.0f, 0.0f, 0.0f);
			}

			Colour F = ShadingHelper::fresnelConductor(cosTheta, eta, k);
			Colour albedoSample = albedo->sample(shadingData.tu, shadingData.tv);
			reflectedColour = (F * albedoSample) / cosTheta;
			pdf = 1.0f;
			return shadingData.frame.toWorld(wiLocal);
		}

		float r1 = sampler->next();
		float r2 = sampler->next();

		float phi_h = 2.0f * (float)M_PI * r2;
		float cosTheta_h = sqrtf((1.0f - r1) / ((alpha * alpha * r1) + (1.0f - r1)));
		float sinTheta_h = sqrtf(std::max(0.0f, 1.0f - cosTheta_h * cosTheta_h));

		Vec3 h = Vec3(sinTheta_h * cosf(phi_h), sinTheta_h * sinf(phi_h), cosTheta_h);

		float dot_wo_h = woLocal.dot(h);
		if (dot_wo_h <= 0.0f)
		{
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return Vec3(0.0f, 0.0f, 0.0f);
		}

		Vec3 wiLocal = (h * (2.0f * dot_wo_h)) - woLocal;
		if (wiLocal.z <= 0.0f)
		{
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return Vec3(0.0f, 0.0f, 0.0f);
		}

		float D = ShadingHelper::Dggx(h, alpha);
		pdf = D * h.z / (4.0f * dot_wo_h);

		Vec3 wiWorld = shadingData.frame.toWorld(wiLocal);
		reflectedColour = evaluate(shadingData, wiWorld);
		return wiWorld;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		if (alpha < 0.001f)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}

		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		if (wiLocal.z <= 0.0f || woLocal.z <= 0.0f)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}

		Vec3 hSum = wiLocal + woLocal;
		if (hSum.lengthSq() <= 0.0f)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}
		Vec3 h = hSum.normalize();

		float cosTheta_i = wiLocal.z;
		float cosTheta_o = woLocal.z;
		float cosTheta_h = wiLocal.dot(h);

		Colour F = ShadingHelper::fresnelConductor(cosTheta_h, eta, k);
		float D = ShadingHelper::Dggx(h, alpha);
		float G = ShadingHelper::Gggx(wiLocal, woLocal, alpha);

		Colour albedoSample = albedo->sample(shadingData.tu, shadingData.tv);
		float denom = 4.0f * cosTheta_i * cosTheta_o;
		return (albedoSample * F) * (D * G / denom);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		if (alpha < 0.001f)
		{
			return 0.0f;
		}

		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		Vec3 hSum = wiLocal + woLocal;
		if (hSum.lengthSq() <= 0.0f)
		{
			return 0.0f;
		}
		Vec3 h = hSum.normalize();

		if (h.z <= 0.0f)
		{
			return 0.0f;
		}

		float dot_wo_h = woLocal.dot(h);
		if (dot_wo_h <= 0.0f)
		{
			return 0.0f;
		}

		float D = ShadingHelper::Dggx(h, alpha);
		return D * h.z / (4.0f * dot_wo_h);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class GlassBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	GlassBSDF() = default;
	GlassBSDF(Texture* _albedo, float _intIOR, float _extIOR)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}
	// 7-5 判断入射方向：光线可能从外面射进去，也可能从里面射出来（比如光在玻璃球内部弹射），两种情况折射率的比值是反的 计算菲涅尔F：F是反射的概率，1 - F是透射的概率。用随机数决定走哪条路，这样蒙特卡洛期望等于真实的菲涅尔加权结果
	// sinThetaTSq >= 1：这是全内反射条件，此时折射不存在，必须强制反射
	// reflectedColour乘以(eta_i / eta_t) ^ 2：光线穿过折射率不同的界面时，光锥会被压缩或扩张，能量密度变化，需要这个Jacobian来补偿
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);

		float etaI = (woLocal.z > 0.0f) ? extIOR : intIOR;
		float etaT = (woLocal.z > 0.0f) ? intIOR : extIOR;

		float F = ShadingHelper::fresnelDielectric(fabsf(woLocal.z), intIOR, extIOR);

		Vec3 wiLocal;
		if (sampler->next() < F)
		{
			wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
			pdf = F;
		}
		else
		{
			float eta = etaI / etaT;
			float cosThetaI = fabsf(woLocal.z);
			float sinThetaTSq = eta * eta * (1.0f - cosThetaI * cosThetaI);

			if (sinThetaTSq >= 1.0f)
			{
				wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
				reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
				pdf = F;
			}
			else
			{
				float cosThetaT = sqrtf(1.0f - sinThetaTSq);
				float sign = (woLocal.z > 0.0f) ? -1.0f : 1.0f;
				wiLocal = Vec3(-eta * woLocal.x, -eta * woLocal.y, sign * cosThetaT);
				reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * (eta * eta);
				pdf = 1.0f - F;
			}
		}

		return shadingData.frame.toWorld(wiLocal);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		(void)shadingData;
		(void)wi;
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		(void)shadingData;
		(void)wi;
		return 0.0f;
	}
	bool isPureSpecular()
	{
		return true;
	}
	bool isTwoSided()
	{
		return false;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class DielectricBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	DielectricBSDF() = default;
	DielectricBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Dielectric sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Dielectric evaluation code
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Dielectric PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return false;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class OrenNayarBSDF : public BSDF
{
public:
	Texture* albedo;
	float sigma;
	OrenNayarBSDF() = default;
	OrenNayarBSDF(Texture* _albedo, float _sigma)
	{
		albedo = _albedo;
		sigma = _sigma;
	}
	// 7-10 evaluate：Oren - Nayar把表面建模成V形沟槽的集合，比朗伯漫反射更接近真实粗糙表面（如月球表面、陶土）
		// A项：主要的漫反射贡献，当sigma = 0时A = 1，退化为朗伯
		// B项：描述逆反射现象（光从哪个方向来就往哪个方向反射得最亮），这在粗糙表面很明显
		// cos(phi_i - phi_o)：描述wi和wo是否在同一个方位角平面内，同一平面时B项贡献最大
		// sin(alpha) * tan(beta)：这两个角度描述光路在V形沟槽里的几何关系
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = SamplingDistributions::cosineHemispherePDF(wi);
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);

		if (wiLocal.z <= 0.0f)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}

		float sigmaSq = sigma * sigma;
		float A = 1.0f - sigmaSq / (2.0f * (sigmaSq + 0.33f));
		float B = 0.45f * sigmaSq / (sigmaSq + 0.09f);

		float theta_i = acosf(std::max(-1.0f, std::min(1.0f, wiLocal.z)));
		float theta_o = acosf(std::max(-1.0f, std::min(1.0f, woLocal.z)));
		float alpha = std::max(theta_i, theta_o);
		float beta = std::min(theta_i, theta_o);

		float sinTheta_i = sinf(theta_i);
		float sinTheta_o = sinf(theta_o);

		float cosDeltaPhi = 0.0f;
		if (sinTheta_i > 1e-4f && sinTheta_o > 1e-4f)
		{
			cosDeltaPhi = (wiLocal.x * woLocal.x + wiLocal.y * woLocal.y) / (sinTheta_i * sinTheta_o);
			cosDeltaPhi = std::max(-1.0f, std::min(1.0f, cosDeltaPhi));
		}

		float orenNayarTerm = A + B * std::max(0.0f, cosDeltaPhi) * sinf(alpha) * tanf(beta);
		return (albedo->sample(shadingData.tu, shadingData.tv) / M_PI) * orenNayarTerm;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class PlasticBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	PlasticBSDF() = default;
	PlasticBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}
	float alphaToPhongExponent()
	{
		return (2.0f / SQ(std::max(alpha, 0.001f))) - 2.0f;
	}
	// 7-11塑料 = 漫反射底层 + 高光层，用菲涅尔F来控制两者比例（掠射角时高光更强）。
		//evaluate：
		//漫反射项：和DiffuseBSDF一样，albedo / π
		//高光项：Phong模型，r是镜面反射方向，wi越接近r，高光越亮，指数n控制高光锐度
		//(n + 2) / (2π)：Phong的归一化系数，保证积分 <= 1
		//F加权：符合物理，掠射时高光更强（菲涅尔效应）

		//sample：		
		//先算F决定走哪个分支：物理上菲涅尔F比例的光走镜面，其余走漫反射
		//高光采样：pow(r1, 1 / (n + 1))生成符合cosⁿθ分布的极角，然后在反射方向r建坐标系，把这个方向转过去
		//不同的pdf：两个分支的pdf不同，要分别记录并反映在返回值里
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		float F = ShadingHelper::fresnelDielectric(fabsf(woLocal.z), intIOR, extIOR);

		Vec3 wiLocal;
		if (sampler->next() < F)
		{
			float r1 = sampler->next();
			float r2 = sampler->next();

			float n = alphaToPhongExponent();
			float cosTheta = powf(r1, 1.0f / (n + 1.0f));
			float sinTheta = sqrtf(std::max(0.0f, 1.0f - cosTheta * cosTheta));
			float phi = 2.0f * (float)M_PI * r2;

			Vec3 r = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			Vec3 wiLocalLobe = Vec3(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);

			Frame lobeFrame;
			lobeFrame.fromVector(r);
			wiLocal = lobeFrame.toWorld(wiLocalLobe);

			if (wiLocal.z <= 0.0f)
			{
				reflectedColour = Colour(0.0f, 0.0f, 0.0f);
				pdf = 0.0f;
				return Vec3(0.0f, 0.0f, 0.0f);
			}

			pdf = powf(cosTheta, n) * (n + 1.0f) / (2.0f * (float)M_PI) * F;
		}
		else
		{
			wiLocal = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			pdf = SamplingDistributions::cosineHemispherePDF(wiLocal) * (1.0f - F);
		}

		Vec3 wiWorld = shadingData.frame.toWorld(wiLocal);
		reflectedColour = evaluate(shadingData, wiWorld);
		return wiWorld;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		if (wiLocal.z <= 0.0f)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}

		Colour diffuse = albedo->sample(shadingData.tu, shadingData.tv) / (float)M_PI;

		Vec3 r = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		float cosPhi = std::max(0.0f, r.dot(wiLocal));
		float n = alphaToPhongExponent();
		Colour specular = Colour(1.0f, 1.0f, 1.0f) * powf(cosPhi, n) * ((n + 2.0f) / (2.0f * (float)M_PI));

		float F = ShadingHelper::fresnelDielectric(woLocal.z, intIOR, extIOR);
		return (diffuse * (1.0f - F)) + (specular * F);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		float F = ShadingHelper::fresnelDielectric(fabsf(woLocal.z), intIOR, extIOR);

		float diffusePdf = SamplingDistributions::cosineHemispherePDF(wiLocal);

		Vec3 r = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		float cosPhi = std::max(0.0f, r.dot(wiLocal));
		float n = alphaToPhongExponent();
		float specularPdf = powf(cosPhi, n) * (n + 1.0f) / (2.0f * (float)M_PI);

		return diffusePdf * (1.0f - F) + specularPdf * F;
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class LayeredBSDF : public BSDF
{
public:
	BSDF* base;
	Colour sigmaa;
	float thickness;
	float intIOR;
	float extIOR;
	LayeredBSDF() = default;
	LayeredBSDF(BSDF* _base, Colour _sigmaa, float _thickness, float _intIOR, float _extIOR)
	{
		base = _base;
		sigmaa = _sigmaa;
		thickness = _thickness;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}
	//模拟有涂层的材质，比如车漆（金属底色+透明清漆层）、上了亮光剂的木头。
	//sample：
	//
	//	F决定哪层响应：涂层界面有一定反射率F，这部分光根本进不了底层，直接被涂层表面反射
	//	折射进去再出来：(1 - F) ^ 2是两次透射（进一次出一次）的乘积，表示能穿透涂层的总比例
	//	吸收项exp(...)：光在涂层内传播会被吸收，厚度thickness越大、吸收系数sigmaa越大，颜色越暗。分别对wi和wo方向的路径长度算吸收（因为入射和出射的路径不同）
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		float F = ShadingHelper::fresnelDielectric(fabsf(woLocal.z), intIOR, extIOR);

		if (sampler->next() < F)
		{
			Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			float cosTheta = std::max(1e-4f, fabsf(wiLocal.z));
			reflectedColour = Colour(F, F, F) / cosTheta;
			pdf = F;
			return shadingData.frame.toWorld(wiLocal);
		}
		else
		{
			auto refractLocal = [](const Vec3& v, float etaI, float etaT, Vec3& vt)->bool
				{
					float eta = etaI / etaT;
					float cosThetaI = fabsf(v.z);
					float sinThetaTSq = eta * eta * (1.0f - cosThetaI * cosThetaI);
					if (sinThetaTSq >= 1.0f)
					{
						return false;
					}

					float cosThetaT = sqrtf(std::max(0.0f, 1.0f - sinThetaTSq));
					float sign = (v.z >= 0.0f) ? 1.0f : -1.0f;
					vt = Vec3(-eta * v.x, -eta * v.y, -sign * cosThetaT);
					return true;
				};

			Vec3 woInside;
			if (!refractLocal(woLocal, extIOR, intIOR, woInside))
			{
				Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
				float cosTheta = std::max(1e-4f, fabsf(wiLocal.z));
				reflectedColour = Colour(F, F, F) / cosTheta;
				pdf = F;
				return shadingData.frame.toWorld(wiLocal);
			}

			ShadingData insideData = shadingData;
			insideData.wo = shadingData.frame.toWorld(woInside);

			Colour baseReflectedColour;
			float dummyBasePdf = 0.0f;
			Vec3 wiInsideWorld = base->sample(insideData, sampler, baseReflectedColour, dummyBasePdf);
			Vec3 wiInside = shadingData.frame.toLocal(wiInsideWorld);

			Vec3 wiLocal;
			if (!refractLocal(wiInside, intIOR, extIOR, wiLocal))
			{
				reflectedColour = Colour(0.0f, 0.0f, 0.0f);
				pdf = 0.0f;
				return Vec3(0.0f, 0.0f, 0.0f);
			}

			float absWiInsideZ = std::max(1e-4f, fabsf(wiInside.z));
			float absWoInsideZ = std::max(1e-4f, fabsf(woInside.z));

			Colour attIn(
				expf(-sigmaa.r * thickness / absWoInsideZ),
				expf(-sigmaa.g * thickness / absWoInsideZ),
				expf(-sigmaa.b * thickness / absWoInsideZ));

			Colour attOut(
				expf(-sigmaa.r * thickness / absWiInsideZ),
				expf(-sigmaa.g * thickness / absWiInsideZ),
				expf(-sigmaa.b * thickness / absWiInsideZ));

			Colour att = attIn * attOut;

			float oneMinusF = 1.0f - F;
			reflectedColour = (baseReflectedColour * att) * (oneMinusF * oneMinusF);

			pdf = oneMinusF * base->PDF(insideData, wiInsideWorld);

			return shadingData.frame.toWorld(wiLocal);
		}
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		float F = ShadingHelper::fresnelDielectric(fabsf(woLocal.z), intIOR, extIOR);

		auto refractLocal = [](const Vec3& v, float etaI, float etaT, Vec3& vt)->bool
			{
				float eta = etaI / etaT;
				float cosThetaI = fabsf(v.z);
				float sinThetaTSq = eta * eta * (1.0f - cosThetaI * cosThetaI);
				if (sinThetaTSq >= 1.0f)
				{
					return false;
				}

				float cosThetaT = sqrtf(std::max(0.0f, 1.0f - sinThetaTSq));
				float sign = (v.z >= 0.0f) ? 1.0f : -1.0f;
				vt = Vec3(-eta * v.x, -eta * v.y, -sign * cosThetaT);
				return true;
			};

		Vec3 woInside;
		Vec3 wiInside;
		if (!refractLocal(woLocal, extIOR, intIOR, woInside) || !refractLocal(wiLocal, extIOR, intIOR, wiInside))
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}

		ShadingData insideData = shadingData;
		insideData.wo = shadingData.frame.toWorld(woInside);
		Vec3 wiInsideWorld = shadingData.frame.toWorld(wiInside);

		Colour baseEval = base->evaluate(insideData, wiInsideWorld);

		float absWiInsideZ = std::max(1e-4f, fabsf(wiInside.z));
		float absWoInsideZ = std::max(1e-4f, fabsf(woInside.z));

		Colour attIn(
			expf(-sigmaa.r * thickness / absWoInsideZ),
			expf(-sigmaa.g * thickness / absWoInsideZ),
			expf(-sigmaa.b * thickness / absWoInsideZ));

		Colour attOut(
			expf(-sigmaa.r * thickness / absWiInsideZ),
			expf(-sigmaa.g * thickness / absWiInsideZ),
			expf(-sigmaa.b * thickness / absWiInsideZ));

		Colour att = attIn * attOut;

		float oneMinusF = 1.0f - F;
		Colour baseTerm = (baseEval * att) * (oneMinusF * oneMinusF);

		return baseTerm;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		float F = ShadingHelper::fresnelDielectric(fabsf(woLocal.z), intIOR, extIOR);

		auto refractLocal = [](const Vec3& v, float etaI, float etaT, Vec3& vt)->bool
			{
				float eta = etaI / etaT;
				float cosThetaI = fabsf(v.z);
				float sinThetaTSq = eta * eta * (1.0f - cosThetaI * cosThetaI);
				if (sinThetaTSq >= 1.0f)
				{
					return false;
				}

				float cosThetaT = sqrtf(std::max(0.0f, 1.0f - sinThetaTSq));
				float sign = (v.z >= 0.0f) ? 1.0f : -1.0f;
				vt = Vec3(-eta * v.x, -eta * v.y, -sign * cosThetaT);
				return true;
			};

		Vec3 woInside;
		Vec3 wiInside;
		if (!refractLocal(woLocal, extIOR, intIOR, woInside) || !refractLocal(wiLocal, extIOR, intIOR, wiInside))
		{
			return 0.0f;
		}

		ShadingData insideData = shadingData;
		insideData.wo = shadingData.frame.toWorld(woInside);
		Vec3 wiInsideWorld = shadingData.frame.toWorld(wiInside);

		return (1.0f - F) * base->PDF(insideData, wiInsideWorld);
	}
	bool isPureSpecular()
	{
		return base->isPureSpecular();
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return base->mask(shadingData);
	}
};
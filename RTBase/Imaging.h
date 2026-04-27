#pragma once

#include "Core.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define __STDC_LIB_EXT1__
#include "stb_image_write.h"

// Stop warnings about buffer overruns if size is zero. Size should never be zero and if it is the code handles it.
#pragma warning( disable : 6386)

constexpr float texelScale = 1.0f / 255.0f;

class Texture
{
public:
	Colour* texels;
	float* alpha;
	int width;
	int height;
	int channels;
	void loadDefault()
	{
		width = 1;
		height = 1;
		channels = 3;
		texels = new Colour[1];
		texels[0] = Colour(1.0f, 1.0f, 1.0f);
	}
	void load(std::string filename)
	{
		alpha = NULL;
		if (filename.find(".hdr") != std::string::npos)
		{
			float* textureData = stbi_loadf(filename.c_str(), &width, &height, &channels, 0);
			if (width == 0 || height == 0)
			{
				loadDefault();
				return;
			}
			texels = new Colour[width * height];
			for (int i = 0; i < (width * height); i++)
			{
				texels[i] = Colour(textureData[i * channels], textureData[(i * channels) + 1], textureData[(i * channels) + 2]);
			}
			stbi_image_free(textureData);
			return;
		}
		unsigned char* textureData = stbi_load(filename.c_str(), &width, &height, &channels, 0);
		if (width == 0 || height == 0)
		{
			loadDefault();
			return;
		}
		texels = new Colour[width * height];
		for (int i = 0; i < (width * height); i++)
		{
			texels[i] = Colour(textureData[i * channels] / 255.0f, textureData[(i * channels) + 1] / 255.0f, textureData[(i * channels) + 2] / 255.0f);
		}
		if (channels == 4)
		{
			alpha = new float[width * height];
			for (int i = 0; i < (width * height); i++)
			{
				alpha[i] = textureData[(i * channels) + 3] / 255.0f;
			}
		}
		stbi_image_free(textureData);
	}
	Colour sample(const float tu, const float tv) const
	{
		Colour tex;
		float u = std::max(0.0f, fabsf(tu)) * width;
		float v = std::max(0.0f, fabsf(tv)) * height;
		int x = (int)floorf(u);
		int y = (int)floorf(v);
		float frac_u = u - x;
		float frac_v = v - y;
		float w0 = (1.0f - frac_u) * (1.0f - frac_v);
		float w1 = frac_u * (1.0f - frac_v);
		float w2 = (1.0f - frac_u) * frac_v;
		float w3 = frac_u * frac_v;
		x = x % width;
		y = y % height;
		Colour s[4];
		s[0] = texels[y * width + x];
		s[1] = texels[y * width + ((x + 1) % width)];
		s[2] = texels[((y + 1) % height) * width + x];
		s[3] = texels[((y + 1) % height) * width + ((x + 1) % width)];
		tex = (s[0] * w0) + (s[1] * w1) + (s[2] * w2) + (s[3] * w3);
		return tex;
	}
	float sampleAlpha(const float tu, const float tv) const
	{
		if (alpha == NULL)
		{
			return 1.0f;
		}
		float tex;
		float u = std::max(0.0f, fabsf(tu)) * width;
		float v = std::max(0.0f, fabsf(tv)) * height;
		int x = (int)floorf(u);
		int y = (int)floorf(v);
		float frac_u = u - x;
		float frac_v = v - y;
		float w0 = (1.0f - frac_u) * (1.0f - frac_v);
		float w1 = frac_u * (1.0f - frac_v);
		float w2 = (1.0f - frac_u) * frac_v;
		float w3 = frac_u * frac_v;
		x = x % width;
		y = y % height;
		float s[4];
		s[0] = alpha[y * width + x];
		s[1] = alpha[y * width + ((x + 1) % width)];
		s[2] = alpha[((y + 1) % height) * width + x];
		s[3] = alpha[((y + 1) % height) * width + ((x + 1) % width)];
		tex = (s[0] * w0) + (s[1] * w1) + (s[2] * w2) + (s[3] * w3);
		return tex;
	}
	~Texture()
	{
		delete[] texels;
		if (alpha != NULL)
		{
			delete alpha;
		}
	}
};

class ImageFilter
{
public:
	virtual float filter(const float x, const float y) const = 0;
	virtual int size() const = 0;
};

class BoxFilter : public ImageFilter
{
public:
	float filter(float x, float y) const
	{
		if (fabsf(x) < 0.5f && fabs(y) < 0.5f)
		{
			return 1.0f;
		}
		return 0;
	}
	int size() const
	{
		return 0;
	}
};
// 2-1  当一束光线打中物体并算出了颜色后，我们不能生硬地只把屏幕上对应的那一个方块像素涂上颜色，这样画出来的边缘会有难看的锯齿（马赛克）。**滤波器（Filter）**的作用就像是“晕染”，它把颜色稍微向周围像素散开一点，让画面更平滑代码里已经给了一个最基础的、硬邦邦的“盒子滤波器（Box Filter）”你需要实现更好的。
class GaussianFilter : public ImageFilter
{
public:
	float radius;
	float alpha;
	float expRadius;
	GaussianFilter(float _radius = 2.0f, float _alpha = 2.0f)
	{
		radius = _radius;
		alpha = _alpha;
		expRadius = expf(-alpha * radius * radius);
	}
	float gaussian1D(const float d) const
	{
		float g = expf(-alpha * d * d) - expRadius; // G(d,r)=e^(-a d^2)-e^(-a r^2)
		return (g > 0.0f) ? g : 0.0f;
	}
	float filter(const float x, const float y) const
	{
		if (fabsf(x) > radius || fabsf(y) > radius)
		{
			return 0.0f;
		}
		return gaussian1D(x) * gaussian1D(y);
	}
	int size() const
	{
		return (int)ceilf(radius);
	}
};
// 2-1
class MitchellNetravaliFilter : public ImageFilter
{
public:
	float B;
	float C;
	float radius;
	MitchellNetravaliFilter(float _B = 1.0f / 3.0f, float _C = 1.0f / 3.0f, float _radius = 2.0f)
	{
		B = _B;
		C = _C;
		radius = _radius;
	}
	float mitchell1D(const float x) const
	{
		// Normalize to canonical Mitchell support [-2, 2]
		float t = fabsf(x) * (2.0f / radius);

		if (t < 1.0f)
		{
			return ((12.0f - 9.0f * B - 6.0f * C) * t * t * t +
				(-18.0f + 12.0f * B + 6.0f * C) * t * t +
				(6.0f - 2.0f * B)) / 6.0f;
		}
		if (t < 2.0f)
		{
			return ((-B - 6.0f * C) * t * t * t +
				(6.0f * B + 30.0f * C) * t * t +
				(-12.0f * B - 48.0f * C) * t +
				(8.0f * B + 24.0f * C)) / 6.0f;
		}
		return 0.0f;
	}
	float filter(const float x, const float y) const
	{
		if (fabsf(x) > radius || fabsf(y) > radius)
		{
			return 0.0f;
		}
		return mitchell1D(x) * mitchell1D(y);
	}
	int size() const
	{
		return (int)ceilf(radius);
	}
};

class Film
{
public:
	Colour* film;
	unsigned int width;
	unsigned int height;
	int SPP;
	ImageFilter* filter;
	// 2-2  在现实里，相机的底片（Film）是用来不断接收光线的。这部分要求你在代码里模拟一张“虚拟胶片（Film Class）”.当你用刚才做好的滤波器算出了光线该怎么“晕染”后，你需要把这些带权重的颜色真正地记录、累加到这张虚拟胶片上。
	void splat(const float x, const float y, const Colour& L)
	{
		int radius = filter->size();

		int minX = (int)floorf(x) - radius;
		int maxX = (int)floorf(x) + radius;
		int minY = (int)floorf(y) - radius;
		int maxY = (int)floorf(y) + radius;

		for (int py = minY; py <= maxY; py++)
		{
			if (py < 0 || py >= (int)height)
			{
				continue;
			}
			for (int px = minX; px <= maxX; px++)
			{
				if (px < 0 || px >= (int)width)
				{
					continue;
				}

				// Sample position relative to pixel center
				float dx = ((float)px + 0.5f) - x;
				float dy = ((float)py + 0.5f) - y;

				float w = filter->filter(dx, dy);
				if (w <= 0.0f)
				{
					continue;
				}

				unsigned int index = (unsigned int)py * width + (unsigned int)px;
				film[index] = film[index] + (L * w);
			}
		}
	}
	// 2-3 虚拟胶片现在记录的亮度可能高达好几万（比如太阳的高光），但你的普通显示器最高只能显示 255 的纯白。如果你直接把高光输出到屏幕，画面就会白茫茫一片，什么细节都看不见。** 色调映射（Tone Mapping）** 就是一门把几万的亮度“温和地压缩”到 0 - 255 范围内，同时还能保留亮部和暗部细节的技术
	void tonemap(int x, int y, unsigned char& r, unsigned char& g, unsigned char& b, float exposure = 1.0f)
	{
		// Safety
		if (x < 0 || x >= (int)width || y < 0 || y >= (int)height)
		{
			r = g = b = 0;
			return;
		}

		unsigned int idx = (unsigned int)y * width + (unsigned int)x;

		// Convert accumulated film value to average radiance
		float invSPP = (SPP > 0) ? (1.0f / (float)SPP) : 1.0f;
		float inR = film[idx].r * invSPP * exposure;
		float inG = film[idx].g * invSPP * exposure;
		float inB = film[idx].b * invSPP * exposure;

		// Clamp negative values before tonemapping
		inR = std::max(0.0f, inR);
		inG = std::max(0.0f, inG);
		inB = std::max(0.0f, inB);

		// Reinhard global operator (normalized by white point)
		// Lout = (C(Lin) / C(W))^(1/2.2), where C(x)=x/(1+x)
		const float whitePoint = 4.0f;
		const float Cw = whitePoint / (1.0f + whitePoint);

		float outR = (inR / (1.0f + inR)) / Cw;
		float outG = (inG / (1.0f + inG)) / Cw;
		float outB = (inB / (1.0f + inB)) / Cw;

		// Gamma correction
		const float invGamma = 1.0f / 2.2f;
		outR = powf(std::max(0.0f, outR), invGamma);
		outG = powf(std::max(0.0f, outG), invGamma);
		outB = powf(std::max(0.0f, outB), invGamma);

		// Clamp to display range [0,1]
		outR = std::min(1.0f, std::max(0.0f, outR));
		outG = std::min(1.0f, std::max(0.0f, outG));
		outB = std::min(1.0f, std::max(0.0f, outB));

		r = (unsigned char)(outR * 255.0f);
		g = (unsigned char)(outG * 255.0f);
		b = (unsigned char)(outB * 255.0f);
	}
	// Do not change any code below this line
	void init(int _width, int _height, ImageFilter* _filter)
	{
		width = _width;
		height = _height;
		film = new Colour[width * height];
		clear();
		filter = _filter;
	}
	void clear()
	{
		memset(film, 0, width * height * sizeof(Colour));
		SPP = 0;
	}
	void incrementSPP()
	{
		SPP++;
	}
	void save(std::string filename)
	{
		Colour* hdrpixels = new Colour[width * height];
		for (unsigned int i = 0; i < (width * height); i++)
		{
			hdrpixels[i] = film[i] / (float)SPP;
		}
		stbi_write_hdr(filename.c_str(), width, height, 3, (float*)hdrpixels);
		delete[] hdrpixels;
	}
};
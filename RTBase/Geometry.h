#pragma once

#include "Core.h"
#include "Sampling.h"

#define EPSILON 0.001f

class Ray
{
public:
	Vec3 o;
	Vec3 dir;
	Vec3 invDir;
	Ray()
	{
	}
	Ray(Vec3 _o, Vec3 _d)
	{
		init(_o, _d);
	}
	void init(Vec3 _o, Vec3 _d)
	{
		o = _o;
		dir = _d;
		invDir = Vec3(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);
	}
	Vec3 at(const float t) const
	{
		return (o + (dir * t));
	}
};

class Plane
{
public:
	Vec3 n;
	float d;
	void init(Vec3& _n, float _d)
	{
		n = _n;
		d = _d;
	}
	// Add code here：1-2 判断刚才发射出去的光线，有没有撞上虚拟世界里的“无限大平面（Plane）”（比如地板或墙壁）
	bool rayIntersect(Ray& r, float& t)
	{
		float denom = Dot(n, r.dir);
		if (fabsf(denom) < EPSILON)
		{
			return false; // Ray is parallel to the plane
		}

		t = (d - Dot(n, r.o)) / denom;
		if (t <= EPSILON)
		{
			return false; // Intersection is behind the ray origin or too close
		}

		return true;
	}
};



class Triangle
{
public:
	Vertex vertices[3];
	Vec3 e1; // Edge 1
	Vec3 e2; // Edge 2
	Vec3 n; // Geometric Normal
	float area; // Triangle area
	float d; // For ray triangle if needed
	unsigned int materialIndex;
	void init(Vertex v0, Vertex v1, Vertex v2, unsigned int _materialIndex)
	{
		materialIndex = _materialIndex;
		vertices[0] = v0;
		vertices[1] = v1;
		vertices[2] = v2;
		e1 = vertices[2].p - vertices[1].p;
		e2 = vertices[0].p - vertices[2].p;
		n = e1.cross(e2).normalize();
		area = e1.cross(e2).length() * 0.5f;
		d = Dot(n, vertices[0].p);
	}
	Vec3 centre() const
	{
		return (vertices[0].p + vertices[1].p + vertices[2].p) / 3.0f;
	}
	// Add code here 1-4 三角形（Triangle）是复杂3D模型（比如人物、汽车）的基本拼装单元。这一步要判断光线是否打中了某个三角形
	bool rayIntersect(const Ray& r, float& t, float& u, float& v) const
	{
		const Vec3& v0 = vertices[0].p;
		const Vec3& v1 = vertices[1].p;
		const Vec3& v2 = vertices[2].p;

		Vec3 edge1 = v1 - v0;
		Vec3 edge2 = v2 - v0;

		Vec3 pvec = r.dir.cross(edge2);
		float det = Dot(edge1, pvec);

		// Ray parallel to triangle plane
		if (fabsf(det) < EPSILON)
		{
			return false;
		}

		// Optimization: replace divisions by multiplication with inverse determinant
		float invDet = 1.0f / det;

		Vec3 tvec = r.o - v0;
		float uStd = Dot(tvec, pvec) * invDet;
		if (uStd < 0.0f || uStd > 1.0f)
		{
			return false;
		}

		Vec3 qvec = tvec.cross(edge1);
		float vStd = Dot(r.dir, qvec) * invDet;
		if (vStd < 0.0f || (uStd + vStd) > 1.0f)
		{
			return false;
		}

		float tHit = Dot(edge2, qvec) * invDet;
		if (tHit <= EPSILON)
		{
			return false;
		}

		t = tHit;

		// Remap barycentric outputs to match existing pipeline:
		// alpha=u -> vertex0, beta=v -> vertex1, gamma=1-u-v -> vertex2
		// Standard MT gives: v0 weight = 1-uStd-vStd, v1 weight = uStd, v2 weight = vStd
		u = 1.0f - uStd - vStd; // alpha (vertex0)
		v = uStd;               // beta  (vertex1)

		return true;
	}
	void interpolateAttributes(const float alpha, const float beta, const float gamma, Vec3& interpolatedNormal, float& interpolatedU, float& interpolatedV) const
	{
		interpolatedNormal = vertices[0].normal * alpha + vertices[1].normal * beta + vertices[2].normal * gamma;
		interpolatedNormal = interpolatedNormal.normalize();
		interpolatedU = vertices[0].u * alpha + vertices[1].u * beta + vertices[2].u * gamma;
		interpolatedV = vertices[0].v * alpha + vertices[1].v * beta + vertices[2].v * gamma;
	}
	//7-13 均匀地在三角形面积上随机取一个点，用于面光源采样（需要知道光从哪个位置发出来）。sqrt(r1)和(1 - sqrt(r1))：直接用r1和r2作为重心坐标会在三角形一角产生聚集，先对r1开方可以保证均匀分布在三角形上 u、v、w是重心坐标，三者之和为1，用它们对三个顶点加权就得到三角形内的随机点pdf = 1 / area：均匀采样，概率密度是面积的倒数
	Vec3 sample(Sampler* sampler, float& pdf)
	{
		float r1 = sampler->next();
		float r2 = sampler->next();

		float sqrtR1 = sqrtf(r1);
		float u = 1.0f - sqrtR1;
		float v = r2 * sqrtR1;
		float w = 1.0f - u - v;

		Vec3 p = vertices[0].p * u + vertices[1].p * v + vertices[2].p * w;
		pdf = 1.0f / area;

		return p;
	}
	Vec3 gNormal()
	{
		return (n * (Dot(vertices[0].normal, n) > 0 ? 1.0f : -1.0f));
	}
};

class AABB
{
public:
	Vec3 max;
	Vec3 min;
	AABB()
	{
		reset();
	}
	void reset()
	{
		max = Vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		min = Vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	}
	void extend(const Vec3 p)
	{
		max = Max(max, p);
		min = Min(min, p);
	}
	// Add code here 1-5 AABB 全称是“轴对齐包围盒（Axis Aligned Bounding Box）”。如果一个3D模型有一万个三角形，光线每走一步都要和一万个三角形算一次求交，电脑会卡死。AABB 就是给这个模型套一个简单的方盒子。如果光线连方盒子都没打中，就不需要去算里面的三角形了，这是一种极其重要的加速手段
	bool rayAABB(const Ray& r, float& t)
	{
		float tmin = (min.x - r.o.x) * r.invDir.x;
		float tmax = (max.x - r.o.x) * r.invDir.x;
		if (tmin > tmax)
		{
			std::swap(tmin, tmax);
		}

		float tymin = (min.y - r.o.y) * r.invDir.y;
		float tymax = (max.y - r.o.y) * r.invDir.y;
		if (tymin > tymax)
		{
			std::swap(tymin, tymax);
		}

		if ((tmin > tymax) || (tymin > tmax))
		{
			return false;
		}

		if (tymin > tmin)
		{
			tmin = tymin;
		}
		if (tymax < tmax)
		{
			tmax = tymax;
		}

		float tzmin = (min.z - r.o.z) * r.invDir.z;
		float tzmax = (max.z - r.o.z) * r.invDir.z;
		if (tzmin > tzmax)
		{
			std::swap(tzmin, tzmax);
		}

		if ((tmin > tzmax) || (tzmin > tmax))
		{
			return false;
		}

		if (tzmin > tmin)
		{
			tmin = tzmin;
		}
		if (tzmax < tmax)
		{
			tmax = tzmax;
		}

		// Box is behind the ray
		if (tmax <= EPSILON)
		{
			return false;
		}

		// First valid hit distance: entry if outside, exit if inside
		t = (tmin > EPSILON) ? tmin : tmax;
		return true;
	}
	// Add code here 1-5
	bool rayAABB(const Ray& r)
	{
		float t;
		return rayAABB(r, t);
	}
	// Add code here
	float area()
	{
		Vec3 size = max - min;
		return ((size.x * size.y) + (size.y * size.z) + (size.x * size.z)) * 2.0f;
	}
};

class Sphere
{
public:
	Vec3 centre;
	float radius;
	void init(Vec3& _centre, float _radius)
	{
		centre = _centre;
		radius = _radius;
	}
	// Add code here:1-3 判断光线有没有打中3D世界里的完美球体（Sphere）
	bool rayIntersect(Ray& r, float& t)
	{
		Vec3 oc = r.o - centre;

		float a = Dot(r.dir, r.dir);
		float b = 2.0f * Dot(oc, r.dir);
		float c = Dot(oc, oc) - (radius * radius);

		float discriminant = (b * b) - (4.0f * a * c);
		if (discriminant < 0.0f)
		{
			return false; // No real roots => miss
		}

		float inv2a = 1.0f / (2.0f * a);

		// Tangent (one root) or two roots
		if (fabsf(discriminant) < EPSILON)
		{
			float t0 = (-b) * inv2a;
			if (t0 > EPSILON)
			{
				t = t0;
				return true;
			}
			return false;
		}

		float sqrtDisc = sqrtf(discriminant);
		float t0 = (-b - sqrtDisc) * inv2a; // near root
		float t1 = (-b + sqrtDisc) * inv2a; // far root

		// Pick nearest valid positive t.
		// If ray starts inside sphere, t0 is usually negative and t1 is the exit hit.
		if (t0 > EPSILON)
		{
			t = t0;
			return true;
		}
		if (t1 > EPSILON)
		{
			t = t1;
			return true;
		}

		return false; // both behind origin
	}
};

struct IntersectionData
{
	unsigned int ID;
	float t;
	float alpha;
	float beta;
	float gamma;
};

#define MAXNODE_TRIANGLES 8
#define TRAVERSE_COST 1.0f
#define TRIANGLE_COST 2.0f
#define BUILD_BINS 32
// 3-1 你要在代码里把“字典的目录”建起来。BVH 的原理是用大盒子包住小盒子，小盒子再包住三角形。你需要定义这个“盒子节点”长什么样
class BVHNode
{
public:
	AABB bounds;
	BVHNode* r;
	BVHNode* l;

	// Leaf storage in reordered global triangle array:
	// [firstTriangle, firstTriangle + triangleCount)
	unsigned int firstTriangle;
	unsigned int triangleCount;

	BVHNode()
	{
		r = NULL;
		l = NULL;
		firstTriangle = 0;
		triangleCount = 0;
	}

	bool isLeaf() const
	{
		return (l == NULL && r == NULL);
	}

	void build(std::vector<Triangle>& inputTriangles)
	{
		r = NULL;
		l = NULL;
		firstTriangle = 0;
		triangleCount = 0;
		if (inputTriangles.empty())
		{
			bounds.reset();
			return;
		}
		buildRecursive(inputTriangles, 0u, (unsigned int)inputTriangles.size());
	}
	//3-3 树建好了，现在要让相机射出的光线去“查字典”了。你要写一段代码，让光线顺着你刚刚建好的层层盒子去找它到底撞到了哪个三角形
	void traverse(const Ray& ray, const std::vector<Triangle>& triangles, IntersectionData& intersection)
	{
		float tNode;
		if (!bounds.rayAABB(ray, tNode))
		{
			return;
		}
		// 如果这个节点最近可能命中的距离都比当前已知最近交点更远，直接剪枝
		if (tNode > intersection.t)
		{
			return;
		}

		// 叶子节点：和该叶子覆盖的三角形做求交
		if (isLeaf())
		{
			unsigned int end = firstTriangle + triangleCount;
			if (end > (unsigned int)triangles.size())
			{
				end = (unsigned int)triangles.size();
			}
			for (unsigned int i = firstTriangle; i < end; i++)
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
			return;
		}

		// 内部节点：先测左右子盒，并按更近的盒子先遍历
		float tL = FLT_MAX;
		float tR = FLT_MAX;
		bool hitL = (l != NULL) && l->bounds.rayAABB(ray, tL);
		bool hitR = (r != NULL) && r->bounds.rayAABB(ray, tR);

		if (hitL && hitR)
		{
			if (tL < tR)
			{
				l->traverse(ray, triangles, intersection);
				r->traverse(ray, triangles, intersection);
			} else
			{
				r->traverse(ray, triangles, intersection);
				l->traverse(ray, triangles, intersection);
			}
		} else if (hitL)
		{
			l->traverse(ray, triangles, intersection);
		} else if (hitR)
		{
			r->traverse(ray, triangles, intersection);
		}
	}

	IntersectionData traverse(const Ray& ray, const std::vector<Triangle>& triangles)
	{
		IntersectionData intersection;
		intersection.t = FLT_MAX;
		traverse(ray, triangles, intersection);
		return intersection;
	}
	// 3-4 我们在画阴影的时候，需要从物体表面向光源发射一条“阴影射线”，看看中间有没有被别的物体挡住。如果有挡住，那就是阴影。这里的逻辑和上一步求交略有不同，讲义专门把它拿出来作为一个优化项
	bool traverseVisible(const Ray& ray, const std::vector<Triangle>& triangles, const float maxT)
	{
		// No occlusion needed for invalid/very short segment
		if (maxT <= EPSILON)
		{
			return true;
		}

		float tNode;
		if (!bounds.rayAABB(ray, tNode))
		{
			return true; // This node is not hit, so it cannot occlude
		}

		// Leaf: test contained triangles, early out on first blocker
		if (isLeaf())
		{
			unsigned int end = firstTriangle + triangleCount;
			if (end > (unsigned int)triangles.size())
			{
				end = (unsigned int)triangles.size();
			}
			for (unsigned int i = firstTriangle; i < end; i++)
			{
				float t;
				float u;
				float v;
				if (triangles[i].rayIntersect(ray, t, u, v))
				{
					if (t > EPSILON && t < maxT)
					{
						return false; // Blocked
					}
				}
			}
			return true;
		}

		// Internal node: traverse near child first, early out as soon as blocked
		float tL = FLT_MAX;
		float tR = FLT_MAX;
		bool hitL = (l != NULL) && l->bounds.rayAABB(ray, tL);
		bool hitR = (r != NULL) && r->bounds.rayAABB(ray, tR);

		if (hitL && hitR)
		{
			if (tL < tR)
			{
				if (!l->traverseVisible(ray, triangles, maxT))
				{
					return false;
				}
				return r->traverseVisible(ray, triangles, maxT);
			} else
			{
				if (!r->traverseVisible(ray, triangles, maxT))
				{
					return false;
				}
				return l->traverseVisible(ray, triangles, maxT);
			}
		}
		if (hitL)
		{
			return l->traverseVisible(ray, triangles, maxT);
		}
		if (hitR)
		{
			return r->traverseVisible(ray, triangles, maxT);
		}

		return true;
	}

private:
	struct SAHBin
	{
		AABB b;
		unsigned int count;
		SAHBin()
		{
			count = 0;
		}
	};

	static float axisValue(const Vec3& v, const int axis)
	{
		if (axis == 0)
		{
			return v.x;
		}
		if (axis == 1)
		{
			return v.y;
		}
		return v.z;
	}

	static float centroidAxis(const Triangle& t, const int axis)
	{
		Vec3 c = t.centre();
		return axisValue(c, axis);
	}

	static AABB triangleBounds(const Triangle& t)
	{
		AABB b;
		b.extend(t.vertices[0].p);
		b.extend(t.vertices[1].p);
		b.extend(t.vertices[2].p);
		return b;
	}

	void buildRecursive(std::vector<Triangle>& inputTriangles, const unsigned int start, const unsigned int end)
	{
		firstTriangle = start;
		triangleCount = end - start;
		l = NULL;
		r = NULL;

		// Compute node bounds + centroid bounds
		bounds.reset();
		AABB centroidBounds;
		centroidBounds.reset();
		for (unsigned int i = start; i < end; i++)
		{
			AABB tb = triangleBounds(inputTriangles[i]);
			bounds.extend(tb.min);
			bounds.extend(tb.max);
			Vec3 c = inputTriangles[i].centre();
			centroidBounds.extend(c);
		}

		// Leaf stop condition
		if (triangleCount <= MAXNODE_TRIANGLES)
		{
			return;
		}

		const float parentArea = bounds.area();
		if (parentArea <= EPSILON)
		{
			return;
		}

		const float leafCost = (float)triangleCount * TRIANGLE_COST;
		float bestCost = FLT_MAX;
		int bestAxis = -1;
		int bestSplit = -1;

		// Evaluate SAH splits on X/Y/Z
		for (int axis = 0; axis < 3; axis++)
		{
			const float cmin = axisValue(centroidBounds.min, axis);
			const float cmax = axisValue(centroidBounds.max, axis);
			const float extent = cmax - cmin;
			if (extent <= EPSILON)
			{
				continue;
			}

			SAHBin bins[BUILD_BINS];

			// Fill bins
			for (unsigned int i = start; i < end; i++)
			{
				const Triangle& t = inputTriangles[i];
				float c = centroidAxis(t, axis);
				float n = (c - cmin) / extent; // [0,1]
				int b = (int)(n * (float)BUILD_BINS);
				if (b < 0)
				{
					b = 0;
				}
				if (b >= BUILD_BINS)
				{
					b = BUILD_BINS - 1;
				}
				bins[b].count++;
				AABB tb = triangleBounds(t);
				bins[b].b.extend(tb.min);
				bins[b].b.extend(tb.max);
			}

			// Prefix/suffix accumulation
			AABB leftBounds[BUILD_BINS - 1];
			AABB rightBounds[BUILD_BINS - 1];
			unsigned int leftCount[BUILD_BINS - 1] = {};
			unsigned int rightCount[BUILD_BINS - 1] = {};

			AABB runLeft;
			runLeft.reset();
			unsigned int runLeftCount = 0;
			for (int i = 0; i < BUILD_BINS - 1; i++)
			{
				runLeftCount += bins[i].count;
				if (bins[i].count > 0)
				{
					runLeft.extend(bins[i].b.min);
					runLeft.extend(bins[i].b.max);
				}
				leftCount[i] = runLeftCount;
				leftBounds[i] = runLeft;
			}

			AABB runRight;
			runRight.reset();
			unsigned int runRightCount = 0;
			for (int i = BUILD_BINS - 1; i > 0; i--)
			{
				runRightCount += bins[i].count;
				if (bins[i].count > 0)
				{
					runRight.extend(bins[i].b.min);
					runRight.extend(bins[i].b.max);
				}
				const int splitIndex = i - 1;
				rightCount[splitIndex] = runRightCount;
				rightBounds[splitIndex] = runRight;
			}

			// Evaluate split cost
			for (int s = 0; s < BUILD_BINS - 1; s++)
			{
				if (leftCount[s] == 0 || rightCount[s] == 0)
				{
					continue;
				}

				const float leftArea = leftBounds[s].area();
				const float rightArea = rightBounds[s].area();

				const float cost =
					TRAVERSE_COST +
					(leftArea / parentArea) * ((float)leftCount[s] * TRIANGLE_COST) +
					(rightArea / parentArea) * ((float)rightCount[s] * TRIANGLE_COST);

				if (cost < bestCost)
				{
					bestCost = cost;
					bestAxis = axis;
					bestSplit = s;
				}
			}
		}

		// If SAH says no benefit, make leaf
		if (bestAxis < 0 || bestCost >= leafCost)
		{
			return;
		}

		// Partition by selected split plane
		const float cmin = axisValue(centroidBounds.min, bestAxis);
		const float cmax = axisValue(centroidBounds.max, bestAxis);
		const float extent = cmax - cmin;
		if (extent <= EPSILON)
		{
			return;
		}
		const float splitPos = cmin + (extent * ((float)(bestSplit + 1) / (float)BUILD_BINS));

		auto beginIt = inputTriangles.begin() + start;
		auto endIt = inputTriangles.begin() + end;

		auto midIt = std::partition(beginIt, endIt,
			[bestAxis, splitPos](const Triangle& t)
			{
				return BVHNode::centroidAxis(t, bestAxis) < splitPos;
			});

		unsigned int mid = (unsigned int)(midIt - inputTriangles.begin());

		// Fallback for degenerate partition
		if (mid == start || mid == end)
		{
			mid = start + (triangleCount / 2);
			auto midNth = inputTriangles.begin() + mid;
			std::nth_element(beginIt, midNth, endIt,
				[bestAxis](const Triangle& a, const Triangle& b)
				{
					return BVHNode::centroidAxis(a, bestAxis) < BVHNode::centroidAxis(b, bestAxis);
				});
		}

		if (mid == start || mid == end)
		{
			return;
		}

		// Internal node
		l = new BVHNode();
		r = new BVHNode();
		l->buildRecursive(inputTriangles, start, mid);
		r->buildRecursive(inputTriangles, mid, end);

		firstTriangle = 0;
		triangleCount = 0;
	}
};

#include "GEMLoader.h"
#include "Renderer.h"
#include "SceneLoader.h"
#include "Geometry.h"
#define NOMINMAX
#include "GamesEngineeringBase.h"
#include <unordered_map>

void runTests()
{
	std::cout << "Running Plane::rayIntersect tests..." << std::endl;

	Plane plane;
	Vec3 n(0.0f, 1.0f, 0.0f);
	plane.init(n, 0.0f); // y = 0 plane

	// Test 1: hit
	{
		Ray ray(Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f));
		float t = -1.0f;
		bool hit = plane.rayIntersect(ray, t);
		bool pass = hit && fabsf(t - 1.0f) < 1e-4f;
		std::cout << "[Plane Test 1] hit from above: " << (pass ? "PASS" : "FAIL") << " (t=" << t << ")" << std::endl;
	}

	// Test 2: parallel
	{
		Ray ray(Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f));
		float t = -1.0f;
		bool hit = plane.rayIntersect(ray, t);
		std::cout << "[Plane Test 2] parallel ray: " << (!hit ? "PASS" : "FAIL") << std::endl;
	}

	// Test 3: intersection behind origin
	{
		Ray ray(Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f));
		float t = -1.0f;
		bool hit = plane.rayIntersect(ray, t);
		std::cout << "[Plane Test 3] behind origin: " << (!hit ? "PASS" : "FAIL") << " (t=" << t << ")" << std::endl;
	}

	std::cout << "Running Sphere::rayIntersect tests..." << std::endl;

	Sphere sphere;
	Vec3 c(0.0f, 0.0f, 0.0f);
	sphere.init(c, 1.0f);

	// Test 1: two intersections (outside -> hit near root)
	{
		Ray ray(Vec3(0.0f, 0.0f, -3.0f), Vec3(0.0f, 0.0f, 1.0f));
		float t = -1.0f;
		bool hit = sphere.rayIntersect(ray, t);
		bool pass = hit && fabsf(t - 2.0f) < 1e-4f;
		std::cout << "[Sphere Test 1] outside hit (2 roots): " << (pass ? "PASS" : "FAIL") << " (t=" << t << ")" << std::endl;
	}

	// Test 2: tangent (discriminant == 0)
	{
		Ray ray(Vec3(1.0f, 0.0f, -3.0f), Vec3(0.0f, 0.0f, 1.0f));
		float t = -1.0f;
		bool hit = sphere.rayIntersect(ray, t);
		bool pass = hit && fabsf(t - 3.0f) < 1e-4f;
		std::cout << "[Sphere Test 2] tangent hit (1 root): " << (pass ? "PASS" : "FAIL") << " (t=" << t << ")" << std::endl;
	}

	// Test 3: miss (discriminant < 0)
	{
		Ray ray(Vec3(2.0f, 0.0f, -3.0f), Vec3(0.0f, 0.0f, 1.0f));
		float t = -1.0f;
		bool hit = sphere.rayIntersect(ray, t);
		std::cout << "[Sphere Test 3] miss: " << (!hit ? "PASS" : "FAIL") << std::endl;
	}

	// Test 4: origin inside sphere
	{
		Ray ray(Vec3(0.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f));
		float t = -1.0f;
		bool hit = sphere.rayIntersect(ray, t);
		bool pass = hit && fabsf(t - 1.0f) < 1e-4f;
		std::cout << "[Sphere Test 4] inside sphere exit hit: " << (pass ? "PASS" : "FAIL") << " (t=" << t << ")" << std::endl;
	}

	std::cout << "Running AABB::rayAABB tests..." << std::endl;

	AABB box;
	box.min = Vec3(-1.0f, -1.0f, -1.0f);
	box.max = Vec3(1.0f, 1.0f, 1.0f);

	// Test 1: outside -> hit (entry)
	{
		Ray ray(Vec3(0.0f, 0.0f, -3.0f), Vec3(0.0f, 0.0f, 1.0f));
		float t = -1.0f;
		bool hit = box.rayAABB(ray, t);
		bool pass = hit && fabsf(t - 2.0f) < 1e-4f;
		std::cout << "[AABB Test 1] outside entry hit: " << (pass ? "PASS" : "FAIL") << " (t=" << t << ")" << std::endl;
	}

	// Test 2: miss
	{
		Ray ray(Vec3(3.0f, 3.0f, -3.0f), Vec3(0.0f, 0.0f, 1.0f));
		bool hit = box.rayAABB(ray);
		std::cout << "[AABB Test 2] miss: " << (!hit ? "PASS" : "FAIL") << std::endl;
	}

	// Test 3: origin inside box -> hit (exit)
	{
		Ray ray(Vec3(0.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f));
		float t = -1.0f;
		bool hit = box.rayAABB(ray, t);
		bool pass = hit && fabsf(t - 1.0f) < 1e-4f;
		std::cout << "[AABB Test 3] inside exit hit: " << (pass ? "PASS" : "FAIL") << " (t=" << t << ")" << std::endl;
	}

	// Test 4: area
	{
		float a = box.area();
		bool pass = fabsf(a - 24.0f) < 1e-4f;
		std::cout << "[AABB Test 4] area: " << (pass ? "PASS" : "FAIL") << " (area=" << a << ")" << std::endl;
	}
}

int main(int argc, char *argv[])
{
	// Add call to tests if required
	runTests();
	
	// Initialize default parameters
	//std::string sceneName = "cornell-box";
	std::string sceneName = "MaterialsScene";
	//std::string sceneName = "bedroom";
	std::string filename = "GI.hdr";
	unsigned int SPP = 8192;

	if (argc > 1)
	{
		std::unordered_map<std::string, std::string> args;
		for (int i = 1; i < argc; ++i)
		{
			std::string arg = argv[i];
			if (!arg.empty() && arg[0] == '-')
			{
				std::string argName = arg;
				if (i + 1 < argc)
				{
					std::string argValue = argv[++i];
					args[argName] = argValue;
				} else
				{
					std::cerr << "Error: Missing value for argument '" << arg << "'\n";
				}
			} else
			{
				std::cerr << "Warning: Ignoring unexpected argument '" << arg << "'\n";
			}
		}
		for (const auto& pair : args)
		{
			if (pair.first == "-scene")
			{
				sceneName = pair.second;
			}
			if (pair.first == "-outputFilename")
			{
				filename = pair.second;
			}
			if (pair.first == "-SPP")
			{
				SPP = stoi(pair.second);
			}
		}
	}
	Scene* scene = loadScene(sceneName);
	// 6-2
	scene->camera.setThinLens(0.001f, (viewcamera.to - viewcamera.from).length());
	GamesEngineeringBase::Window canvas;
	canvas.create((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, "Tracer", false);
	RayTracer rt;
	rt.init(scene, &canvas);
	bool running = true;
	GamesEngineeringBase::Timer timer;
	while (running)
	{
		canvas.checkInput();
		canvas.clear();
		if (canvas.keyPressed(VK_ESCAPE))
		{
			break;
		}
		if (canvas.keyPressed('W'))
		{
			viewcamera.forward();
			rt.clear();
		}
		if (canvas.keyPressed('S'))
		{
			viewcamera.back();
			rt.clear();
		}
		if (canvas.keyPressed('A'))
		{
			viewcamera.left();
			rt.clear();
		}
		if (canvas.keyPressed('D'))
		{
			viewcamera.right();
			rt.clear();
		}
		if (canvas.keyPressed('E'))
		{
			viewcamera.flyUp();
			rt.clear();
		}
		if (canvas.keyPressed('Q'))
		{
			viewcamera.flyDown();
			rt.clear();
		}
		// Time how long a render call takes
		timer.reset();
		rt.render();
		float t = timer.dt();
		// Write
		std::cout << t << std::endl;
		if (canvas.keyPressed('P'))
		{
			rt.saveHDR(filename);
		}
		if (canvas.keyPressed('L'))
		{
			size_t pos = filename.find_last_of('.');
			std::string ldrFilename = filename.substr(0, pos) + ".png";
			rt.savePNG(ldrFilename);
		}
		if (SPP == rt.getSPP())
		{
			rt.saveHDR(filename);
			break;
		}
		canvas.present();
	}
	return 0;
}
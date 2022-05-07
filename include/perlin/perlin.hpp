#pragma once
#include <glm/glm.hpp>
#include <limits>
#include <perlin/analyticPerlin.hpp>

static PerlinNoise analyticPerlin;

float perlin2d(float x, float y, float freq = 0.007f) {
	glm::vec3 p(x * freq, 0, y * freq);
	glm::vec3 d;
	analyticPerlin.eval(p, d);
	return p.y;
}

glm::vec3 perlin2d_normal(float x, float y, float scl, float freq = 0.007f) {
	glm::vec3 p(x * freq, 0, y * freq);
	glm::vec3 d;
	analyticPerlin.eval(p, d);

	return normalize(glm::vec3(-d.x * scl, 1, -d.z * scl));
}
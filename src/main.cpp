
#include <gfx/gfx.hpp>
#include <glfw-helper/glfw.hpp>
#include <perlin/perlin.hpp>

struct app {

	struct Vertex {
		glm::vec4 pos;
		glm::vec4 color;
		glm::vec3 normal;
	};

	int windowWidth = 1000;
	int windowHeight = 600;

	enum BIOME_TYPE { NORMAL, RANDOM };

	BIOME_TYPE biomeType = NORMAL;

	bool firstmouse = true;
	bool randomBiomeFlag = false;
	bool imguiWindowHovered = false;

	glfw::mouse mouse;
	glfw::keyboard keyboard;

	gfx::camera camera;

	gfx::shader basic_shader;
	gfx::uniform<glm::mat4> pvm;
	gfx::uniform<glm::mat4> nm;
	gfx::uniform<glm::vec3> eyePos;

	gfx::buffer<Vertex> terrain;
	std::vector<Vertex> vertices;
	std::vector<unsigned int> indices;
	std::vector<Vertex> modifyVertices;
	std::vector<unsigned int> modifyIndices;

	float scl = 100.0f;
	int resolution = 100;
	int width = 100;
	int height = 100;
	float perlinFreq = 0.01f;

	~app() noexcept { gfx::destroy(terrain); }

	void init() {

		glEnable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);

		basic_shader = gfx::shader::create();
		{
			basic_shader.attach(GL_VERTEX_SHADER, R"(
			#version 330
			
			layout(location = 0) in vec4 pos;
			layout(location = 1) in vec4 col;
			layout(location = 2) in vec3 nor;

			out vec4 color;
			out vec3 worldPos;
			out vec3 normal;

			uniform mat4 pvm;
			uniform mat4 nm;
			uniform vec3 eyePos;

			void main() {
				color = col;
				vec4 position = pvm * pos;
				gl_Position = position;
				worldPos = position.xyz / position.w;
				normal = normalize((nm * vec4(nor, 1.0f)).xyz);
			}

		)");
		}

		{
			basic_shader.attach(GL_FRAGMENT_SHADER, R"(
			#version 330
			
			const vec3 lightDirection = vec3(1.0f, 1.0f, 1.0f);
			const vec3 lightIntensity = vec3(1.0f, 1.0f, 1.0f);

			const vec3 ambientColor = vec3(0.0f, 0.0f, 0.0f);

			const vec3 specColor = vec3(1.0f, 1.0f, 1.0f);
			const float shininess = 16.0f;

			uniform mat4 pvm;
			uniform mat4 nm;
			uniform vec3 eyePos;

			in vec4 color;
			in vec3 worldPos;
			in vec3 normal;

			out vec4 fragColor;

			void main() {
				vec3 lig = normalize(lightDirection);
				vec3 nor = normalize(normal);
				vec3 viewDir = normalize(eyePos - worldPos);
				float lambertian = max(dot(lig, nor), 0.0f);

				vec3 sumColor = vec3(0.0f, 0.0f, 0.0f);

				if (lambertian > 0.0f) {

					vec3 halfDir = normalize(viewDir + lig);
					float cosa = max(dot(halfDir, nor), 0.0f);
					float specular = pow(cosa, shininess);

					sumColor += lightIntensity * (color.rgb * lambertian + specColor * specular);

				} 
				sumColor += ambientColor;
				fragColor = vec4(sumColor, 1.0f);

			}
		)");
		}
		basic_shader.link("fragColor");

		pvm = gfx::uniform<glm::mat4>(basic_shader, "pvm");
		nm = gfx::uniform<glm::mat4>(basic_shader, "nm");
		eyePos = gfx::uniform<glm::vec3>(basic_shader, "eyePos");
		//  pvm = camera.view();

		gfx::layout<Vertex> basic_layout = {
			{0, &Vertex::pos},
			{1, &Vertex::color},
			{2, &Vertex::normal},
		};

		camera = gfx::camera(glm::vec4(1.0f, 70.0f, 3.0f, 1.0f));
		camera.lookat(glm::vec3(width / 2, 1.0f, height / 2));

		generateTerrain(width, height);

		gfx::set_layout(terrain, basic_layout);
	}

	void generateTerrain(int width, int height) {

		vertices.clear();
		indices.clear();

		auto dx = (float)width / (resolution - 1);
		auto dz = (float)height / (resolution - 1);

		auto dw = 1.0f / (resolution - 1.0f) * width;
		auto dh = 1.0f / (resolution - 1.0f) * height;

		for (int i = 0; i < resolution - 1; ++i) {
			for (int j = 0; j < resolution - 1; ++j) {

				const float p_i0_j0 = perlin2d((i + 0) * dw, (j + 0) * dh, perlinFreq);
				const float p_i1_j0 = perlin2d((i + 1) * dw, (j + 0) * dh, perlinFreq);
				const float p_i0_j1 = perlin2d((i + 0) * dw, (j + 1) * dh, perlinFreq);
				const float p_i1_j1 = perlin2d((i + 1) * dw, (j + 1) * dh, perlinFreq);

				glm::vec4 pos1(i * dx, p_i0_j0 * scl, j * dz, 1.0f);
				glm::vec4 pos2((i + 1) * dx, p_i1_j0 * scl, j * dz, 1.0f);
				glm::vec4 pos3(i * dx, p_i0_j1 * scl, (j + 1) * dz, 1.0f);
				glm::vec4 pos4((i + 1) * dx, p_i1_j1 * scl, (j + 1) * dz, 1.0f);

				glm::vec3 norm1 = perlin2d_normal((i + 0) * dw, (j + 0) * dh, scl, perlinFreq);
				glm::vec3 norm2 = perlin2d_normal((i + 1) * dw, (j + 0) * dh, scl, perlinFreq);
				glm::vec3 norm3 = perlin2d_normal((i + 0) * dw, (j + 1) * dh, scl, perlinFreq);
				glm::vec3 norm4 = perlin2d_normal((i + 1) * dw, (j + 1) * dh, scl, perlinFreq);

				glm::vec3 norm1 = perlin2d_normal((i + 0) * dw, (j + 0) * dh, scl, perlinFreq);
				glm::vec3 norm2 = perlin2d_normal((i + 1) * dw, (j + 0) * dh, scl, perlinFreq);
				glm::vec3 norm3 = perlin2d_normal((i + 0) * dw, (j + 1) * dh, scl, perlinFreq);
				glm::vec3 norm4 = perlin2d_normal((i + 1) * dw, (j + 1) * dh, scl, perlinFreq);

				// glm::vec4 cartoonColor1 = {(p_i0_j0 + p_i1_j0 + p_i0_j1) / 3, 0.5f, 0.0f, 1.0f};
				// glm::vec4 cartoonColor2 = {(p_i1_j0 + p_i0_j1 + p_i1_j1) / 3, 0.5f, 0.0f, 1.0f};
				// glm::vec4 randomColor1 = gfx::randVec4();
				// glm::vec4 randomColor2 = gfx::randVec4();
				// randomColor1.a = randomColor2.a = 1.0f;

				vertices.push_back({pos1, heightmapColor(pos1, biomeType)});
				vertices.push_back({pos2, heightmapColor(pos2, biomeType)});
				vertices.push_back({pos3, heightmapColor(pos3, biomeType)});

				/*vertices.push_back({glm::vec4{i * dx, p_i0_j0 * scl, j * dz, 1.0f}, randomColor1});
				vertices.push_back({pos2, randomColor1});
				vertices.push_back({pos3, randomColor1});*/

				vertices.push_back({pos4, heightmapColor(pos4, biomeType)});
				vertices.push_back({pos3, heightmapColor(pos3, biomeType)});
				vertices.push_back({pos2, heightmapColor(pos2, biomeType)});

				/*vertices.push_back({glm::vec4{(i + 1) * dx, p_i1_j1 * scl, (j + 1) * dz, 1.0f}, randomColor2});
				vertices.push_back({pos3, randomColor2});
				vertices.push_back({pos2, randomColor2});*/

				indices.push_back(std::size(indices));
				indices.push_back(std::size(indices));
				indices.push_back(std::size(indices));

				indices.push_back(std::size(indices));
				indices.push_back(std::size(indices));
				indices.push_back(std::size(indices));
			}
		}

		gfx::upload(terrain, gfx::vertex_buffer, vertices);
		gfx::upload(terrain, gfx::index_buffer, indices);
	}

	glm::vec4 heightmapColor(glm::vec4 const &vertex, const BIOME_TYPE biomeType = NORMAL) {
		float colorParam = (perlin2d(vertex.x, vertex.z, perlinFreq) * 50 + vertex.y / scl) / 51;
		switch (biomeType) {
		case NORMAL:
			if (colorParam < 0.55f)
				return glm::vec4(1.0 * colorParam, 0.8f * (1 - colorParam), 0.0f, 1.0f);
			else
				return glm::vec4(1.0, 1.0f, 1.0f, 1.0f);
			break;

		case RANDOM:
			return glm::vec4(gfx::randVec4(0.0f, 1.0f, 0.5f));
			break;

		default:
			return glm::vec4(perlin2d(vertex.x, vertex.z), 0.5f, 0.0f, 1.0f);
			break;
		}
	}

	void key(int key, int scancode, int action, int mods) {
		if (action == GLFW_PRESS) {
			keyboard.keys[key] = true;
		} else if (action == GLFW_RELEASE) {
			keyboard.keys[key] = false;
		}
	}

	void mousebutton(int button, int action, int mode) {
		mouse.buttons[button] = (GLFW_PRESS == action);
		if (mouse.buttons[GLFW_MOUSE_BUTTON_LEFT]) {

			glm::vec3 unprojectNear(0.0f, 0.0f, 0.0f);
			glm::vec3 unprojectFar(0.0f, 0.0f, 0.0f);
			screenToWorldVertex(mouse.lastX, windowHeight - mouse.lastY, unprojectNear, unprojectFar);

			glm::vec3 origin = camera.position;
			glm::vec3 vector = glm::normalize(unprojectNear - unprojectFar);
			Triangle triangle;
			glm::vec3 intersectPoint;
			for (unsigned int i = 0; i < std::size(indices); i += 3) {
				triangle.vertex0 = vertices[i].pos;
				triangle.vertex1 = vertices[i + 1].pos;
				triangle.vertex2 = vertices[i + 2].pos;
				if (rayIntersectsTriangle(origin, vector, triangle, intersectPoint)) {
					std::cout << "Intersection Point [ x: " << intersectPoint[0] << " y: " << intersectPoint[1] << " z: " << intersectPoint[2] << " ]" << std::endl;
					vertices[i].color = glm::vec4(1.0, 0.0, 0.0, 1.0);
					vertices[i + 1].color = glm::vec4(1.0, 0.0, 0.0, 1.0);
					vertices[i + 2].color = glm::vec4(1.0, 0.0, 0.0, 1.0);

					gfx::upload(terrain, gfx::vertex_buffer, vertices);
				}
			}
			// std::cout << rayIntersectsTriangle(origin, vector, &triangle, intersectPoint) << " " << intersectPoint[1];
		}
	}

	void mousemove(double xpos, double ypos) {
		if (mouse.buttons[GLFW_MOUSE_BUTTON_LEFT] && !imguiWindowHovered) {

			if (firstmouse) {
				mouse.lastX = xpos;
				mouse.lastY = ypos;
				firstmouse = false;
			}

			float xoffset = xpos - mouse.lastX;
			float yoffset = mouse.lastY - ypos;
			mouse.lastX = xpos;
			mouse.lastY = ypos;

			xoffset *= mouse.sensitivity;
			yoffset *= mouse.sensitivity;

			camera.yaw -= xoffset;
			camera.pitch -= yoffset;

			if (camera.pitch > 89.0f) camera.pitch = 89.0f;
			if (camera.pitch < -89.0f) camera.pitch = -89.0f;

			glm::vec3 direction;
			direction.x = cos(glm::radians(camera.yaw)) * cos(glm::radians(camera.pitch));
			direction.y = sin(glm::radians(camera.pitch));
			direction.z = sin(glm::radians(camera.yaw)) * cos(glm::radians(camera.pitch));
			camera.front = glm::normalize(direction);

			camera.update_camera_vectors();
		}

		mouse.lastX = xpos;
		mouse.lastY = ypos;
	}

	void screenToWorldVertex(double xpos, double ypos, glm::vec3 &unprojNear, glm::vec3 &unprojFar) {
		unprojNear = glm::unProject(glm::vec3(xpos, ypos, camera.near), camera.view(), camera.projection(), glm::vec4(0, 0, windowWidth, windowHeight));
		unprojFar = glm::unProject(glm::vec3(xpos, ypos, camera.far), camera.view(), camera.projection(), glm::vec4(0, 0, windowWidth, windowHeight));
		// std::cout << "Near: [ " << unprojNear[0] << unprojNear[1] << unprojNear[2] << " ]" << std::endl;
		// std::cout << "Far: [ " << unprojFar[0] << unprojFar[1] << unprojFar[2] << " ]" << std::endl;
		// std::cout << "Size: " << sqrt(std::pow(unprojFar[0] - unprojNear[0], 2) + std::pow(unprojFar[1] - unprojNear[1], 2) + std::pow(unprojFar[2] - unprojNear[2], 2)) << std::endl;
	}

	void resize(int width, int height) {
		windowWidth = width;
		windowHeight = height;
		glViewport(0, 0, width, height);
		camera.aspect = ((float)width) / height;
	}

	void render_gui() {
		auto guiWindowWidth = 300;
		auto guiWindowHeight = 200;
		ImVec2 windowPosition = ImVec2(windowWidth - guiWindowWidth, 0);
		ImVec2 windowSize = ImVec2(guiWindowWidth, guiWindowHeight);

		ImGuiWindowFlags window_flags = 0;
		window_flags |= ImGuiWindowFlags_NoMove;
		window_flags |= ImGuiWindowFlags_NoResize;
		window_flags |= ImGuiWindowFlags_NoCollapse;

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		ImGui::StyleColorsDark();
		ImGui::Begin("Window", NULL, window_flags);
		ImGui::SetWindowPos(windowPosition);
		ImGui::SetWindowSize(windowSize);

		float newScl = scl;
		if (ImGui::SliderFloat("Scale", &newScl, 1.0f, 200.0f)) {
			/*for (auto &vert : vertices)
				vert.pos.y = vert.pos.y / scl * newScl;*/

			auto dx = (float)width / (resolution - 1);
			auto dz = (float)height / (resolution - 1);

			auto dw = 1.0f / (resolution - 1.0f) * width;
			auto dh = 1.0f / (resolution - 1.0f) * height;

			std::size_t k = 0;
			for (int i = 0; i < resolution - 1; ++i) {
				for (int j = 0; j < resolution - 1; ++j) {

					vertices[k].pos.y = vertices[k].pos.y / scl * newScl;
					vertices[k + 1].pos.y = vertices[k + 1].pos.y / scl * newScl;
					vertices[k + 2].pos.y = vertices[k + 2].pos.y / scl * newScl;
					vertices[k + 3].pos.y = vertices[k + 3].pos.y / scl * newScl;
					vertices[k + 4].pos.y = vertices[k + 4].pos.y / scl * newScl;
					vertices[k + 5].pos.y = vertices[k + 5].pos.y / scl * newScl;

					glm::vec3 norm1 = perlin2d_normal((i + 0) * dw, (j + 0) * dh, newScl, perlinFreq);
					glm::vec3 norm2 = perlin2d_normal((i + 1) * dw, (j + 0) * dh, newScl, perlinFreq);
					glm::vec3 norm3 = perlin2d_normal((i + 0) * dw, (j + 1) * dh, newScl, perlinFreq);
					glm::vec3 norm4 = perlin2d_normal((i + 1) * dw, (j + 1) * dh, newScl, perlinFreq);

					vertices[k].normal = norm1;
					vertices[k + 1].normal = norm2;
					vertices[k + 2].normal = norm3;
					vertices[k + 3].normal = norm4;
					vertices[k + 4].normal = norm3;
					vertices[k + 5].normal = norm2;

					k += 6;
				}
			}
			scl = newScl;
			gfx::upload(terrain, gfx::vertex_buffer, vertices);
		}
		// if (ImGui::SliderInt("Resolution", &resolution, 1.0f, 200.0f)) generateTerrain(width, height); TODO
		if (ImGui::SliderFloat("Noise Frequency", &perlinFreq, 0.001f, 0.1f)) generateTerrain(width, height);

		const char *items[] = {"NORMAL", "RANDOM"};
		static int selected_item = 0;
		static int previous_item = 0;

		if (ImGui::Combo("Biome Type", &selected_item, items, IM_ARRAYSIZE(items))) {
			if (selected_item == 0) biomeType = NORMAL;
			if (selected_item == 1) biomeType = RANDOM;
			if (previous_item != selected_item) generateTerrain(width, height);
			previous_item = selected_item;
		}

		imguiWindowHovered = ImGui::IsWindowHovered() || ImGui::IsAnyItemHovered() || ImGui::IsMouseHoveringRect(windowPosition, ImVec2(windowPosition.x + windowSize.x, windowPosition.y + windowSize.y));

		ImGui::End();

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
	}

	void render_terrain() {
		gfx::use(basic_shader);
		pvm = camera.projection() * camera.view();
		nm = glm::mat4(1.0);
		eyePos = camera.position;
		gfx::draw(terrain);
	}

	void render() {
		glClearColor(0, 0, 0, 1);
		glClearDepth(1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		render_terrain();
		render_gui();
	}

	void cameramove(float dt) {
		constexpr auto speed = 10;
		if (keyboard.keys[GLFW_KEY_W]) {
			glm::vec3 m = glm::normalize(glm::vec3{camera.front.x, 0.0f, camera.front.z});

			camera.position.x += m.x * dt * speed;
			camera.position.z += m.z * dt * speed;
		}
		if (keyboard.keys[GLFW_KEY_S]) {
			glm::vec3 m = glm::normalize(glm::vec3{camera.front.x, 0.0f, camera.front.z});

			camera.position.x -= m.x * dt * speed;
			camera.position.z -= m.z * dt * speed;
		}
		if (keyboard.keys[GLFW_KEY_A]) {
			camera.position -= camera.right * dt * speed;
		}
		if (keyboard.keys[GLFW_KEY_D]) {
			camera.position += camera.right * dt * speed;
		}
		if (keyboard.keys[GLFW_KEY_SPACE]) {
			camera.position += camera.worldup * dt * speed;
		}
		if (keyboard.keys[GLFW_KEY_LEFT_SHIFT]) {
			camera.position -= camera.worldup * dt * speed;
		}
	}

	void update(float dt) {
		cameramove(dt);
		render();
	}

	struct Triangle {
		glm::vec3 vertex0;
		glm::vec3 vertex1;
		glm::vec3 vertex2;
	};

	bool rayIntersectsTriangle(glm::vec3 rayOrigin, glm::vec3 rayVector, const Triangle &inTriangle, glm::vec3 &outIntersectionPoint) {
		const float EPSILON = 0.0000001;
		glm::vec3 vertex0 = inTriangle.vertex0;
		glm::vec3 vertex1 = inTriangle.vertex1;
		glm::vec3 vertex2 = inTriangle.vertex2;
		glm::vec3 edge1, edge2, h, s, q;
		float a, f, u, v;
		edge1 = vertex1 - vertex0;
		edge2 = vertex2 - vertex0;
		h = glm::cross(rayVector, edge2);
		a = glm::dot(edge1, h);
		if (a > -EPSILON && a < EPSILON) return false; // This ray is parallel to this triangle.
		f = 1.0 / a;
		s = rayOrigin - vertex0;
		u = f * glm::dot(s, h);
		if (u < 0.0 || u > 1.0) return false;
		q = glm::cross(s, edge1);
		v = f * glm::dot(rayVector, q);
		if (v < 0.0 || u + v > 1.0) return false;
		// At this stage we can compute t to find out where the intersection point is on the line.
		float t = f * glm::dot(edge2, q);
		std::cout << t << std::endl;
		if (t > EPSILON) // ray intersection
		{
			outIntersectionPoint = rayOrigin + rayVector * t;
			return true;
		} else // This means that there is a line intersection but not a ray intersection.
			return false;
	}
};

int main() {
	app app;

	[[maybe_unused]] auto &ctx = glfw::instance();

	glfw::gl::initialize();

	glfw::window wnd(app.windowWidth, app.windowHeight, "WorldGenerator");

	glfw::gl::make_current(wnd);

	glfw::glad::initialize();

	glfw::imgui::initialize(wnd);

	glfw::on(wnd, glfw::resize, [&](GLFWwindow *window, int width, int height) { app.resize(width, height); });
	glfw::on(wnd, glfw::mousemove, [&](glfw::window &window, double xpos, double ypos) { app.mousemove(xpos, ypos); });
	glfw::on(wnd, glfw::mousebutton, [&](glfw::window &window, int button, int action, int mode) { app.mousebutton(button, action, mode); });
	glfw::on(wnd, glfw::key, [&](glfw::window &window, int key, int scancode, int action, int mods) { app.key(key, scancode, action, mods); });

	app.init();

	glfw::run([&](float dt) {
		app.update(dt);

		glfw::swap_buffers(wnd);
		return glfw::is_closed(wnd);
	});

	glfw::imgui::uninitialize();

	return 0;
}
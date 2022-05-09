
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
	bool imguiWindowHovered = false;
	bool brushFlag = false;

	glfw::mouse mouse;
	glfw::keyboard keyboard;

	gfx::camera camera;
	float deltatime;

	gfx::shader basic_shader;
	gfx::shader skybox_shader;

	gfx::uniform<glm::mat4> pvm;
	gfx::uniform<glm::mat4> nm;
	gfx::uniform<glm::vec3> eyePos;

	gfx::buffer<Vertex> terrain;
	gfx::buffer<Vertex> skyBox;
	std::vector<Vertex> vertices;
	std::vector<unsigned int> indices;

	float brushSize = 1.0f;

	float scl = 100.0f;
	const int chunkSize = 500;
	int width = chunkSize;
	int height = chunkSize;
	int resolution = 200;
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
			const float shininess = 8.0f;

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

		skybox_shader = gfx::shader::create();
		{
			skybox_shader.attach(GL_VERTEX_SHADER, R"(
			#version 330
			
			layout(location = 0) in vec4 pos;
			layout(location = 1) in vec4 col;

			out vec4 color;

			void main() {
				color = col;
				gl_Position = pos;
			}

		)");
		}

		{
			skybox_shader.attach(GL_FRAGMENT_SHADER, R"(
			#version 330
			
			in vec4 color;

			out vec4 fragColor;
			
			void main() {
				fragColor = color;
			}
		)");
		}
		skybox_shader.link("fragColor");

		pvm = gfx::uniform<glm::mat4>(basic_shader, "pvm");
		nm = gfx::uniform<glm::mat4>(basic_shader, "nm");
		eyePos = gfx::uniform<glm::vec3>(basic_shader, "eyePos");
		//  pvm = camera.view();

		camera = gfx::camera(glm::vec4(1.0f, 70.0f, 3.0f, 1.0f));
		camera.lookat(glm::vec3(width / 2, 1.0f, height / 2));

		generateTerrain();
		generateSkybox();
	}

	void generateSkybox() {
		std::vector<Vertex> skyboxVertices;
		std::vector<unsigned int> skyboxIndices;

		skyboxVertices.push_back({glm::vec4(1.0f, 1.0f, 0.9999f, 1.0f), glm::vec4(0.52f, 0.82f, 0.99f, 1.0f)});
		skyboxVertices.push_back({glm::vec4(-1.0f, 1.0f, 0.9999f, 1.0f), glm::vec4(0.52f, 0.82f, 0.99f, 1.0f)});
		skyboxVertices.push_back({glm::vec4(1.0f, -1.0f, 0.9999f, 1.0f), glm::vec4(0.52f, 0.82f, 0.99f, 1.0f)});

		skyboxVertices.push_back({glm::vec4(-1.0f, -1.0f, 0.9999f, 1.0f), glm::vec4(0.52f, 0.82f, 0.99f, 1.0f)});
		skyboxVertices.push_back({glm::vec4(1.0f, -1.0f, 0.9999f, 1.0f), glm::vec4(0.52f, 0.82f, 0.99f, 1.0f)});
		skyboxVertices.push_back({glm::vec4(-1.0f, 1.0f, 0.9999f, 1.0f), glm::vec4(0.52f, 0.82f, 0.99f, 1.0f)});

		skyboxIndices.push_back(std::size(skyboxIndices));
		skyboxIndices.push_back(std::size(skyboxIndices));
		skyboxIndices.push_back(std::size(skyboxIndices));

		skyboxIndices.push_back(std::size(skyboxIndices));
		skyboxIndices.push_back(std::size(skyboxIndices));
		skyboxIndices.push_back(std::size(skyboxIndices));

		gfx::upload(skyBox, gfx::vertex_buffer, skyboxVertices);
		gfx::upload(skyBox, gfx::index_buffer, skyboxIndices);
	}

	void generateTerrain() {
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

				vertices.push_back({pos1, heightmapColor(pos1, biomeType), norm1});
				vertices.push_back({pos2, heightmapColor(pos2, biomeType), norm2});
				vertices.push_back({pos3, heightmapColor(pos3, biomeType), norm3});

				vertices.push_back({pos4, heightmapColor(pos4, biomeType), norm4});
				vertices.push_back({pos3, heightmapColor(pos3, biomeType), norm3});
				vertices.push_back({pos2, heightmapColor(pos2, biomeType), norm2});

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

	void terraformingOnClick(short direction) {
		glm::vec3 unprojectNear(0.0f, 0.0f, 0.0f);
		glm::vec3 unprojectFar(0.0f, 0.0f, 0.0f);
		camera.screenToWorldVertex(mouse.lastX, windowHeight - mouse.lastY, unprojectNear, unprojectFar, windowWidth, windowHeight);

		glm::vec3 origin = camera.position;
		glm::vec3 vector = glm::normalize(unprojectNear - unprojectFar);
		gfx::Triangle triangle;
		std::optional<glm::vec3> intersectionPoint;
		float intersectionPointDistance;
		float lastIntersectionPointDistance = std::numeric_limits<float>::max();

		for (unsigned int i = 0; i < std::size(indices); i += 3) {

			triangle.vertex0 = vertices[i].pos;
			triangle.vertex1 = vertices[i + 1].pos;
			triangle.vertex2 = vertices[i + 2].pos;

			if (gfx::rayIntersectsTriangle(origin, vector, triangle, intersectionPointDistance)) {
				if (intersectionPointDistance < lastIntersectionPointDistance) {

					lastIntersectionPointDistance = intersectionPointDistance;
					intersectionPoint = vector * intersectionPointDistance + origin;
				}
			}
		}

		if (intersectionPoint) {
			for (auto &vert : vertices) {
				const float dist = glm::distance(*intersectionPoint, glm::vec3(vert.pos)) / ((width / (resolution - 1) + height / (resolution - 1)) / 2);
				float scale = 1.0f / (1.0f + 1 / brushSize * pow(dist, 2));
				if (scale > 0.1) vert.pos.y += scale * direction;
			}
			gfx::upload(terrain, gfx::vertex_buffer, vertices);
		}
	}

	void rotateCamera(double xpos, double ypos) {
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

	void key(int key, int scancode, int action, int mods) {
		if (action == GLFW_PRESS) {
			keyboard.keys[key] = true;
		} else if (action == GLFW_RELEASE) {
			keyboard.keys[key] = false;
		}
	}

	void mousebutton(int button, int action, int mode) {
		if (action == GLFW_PRESS) {
			mouse.buttons[button] = true;
		} else if (action == GLFW_RELEASE) {
			mouse.buttons[button] = false;
		}
	}

	void mousemove(double xpos, double ypos) {
		if (mouse.buttons[GLFW_MOUSE_BUTTON_LEFT] && !imguiWindowHovered && !brushFlag) {
			rotateCamera(xpos, ypos);
		}

		mouse.lastX = xpos;
		mouse.lastY = ypos;
	}

	void mousewheel(double xoffset, double yoffset) {
		constexpr auto speed = 170;
		camera.position += camera.worldup * deltatime * yoffset * speed;
	}

	void resize(int width, int height) {
		windowWidth = width;
		windowHeight = height;
		glViewport(0, 0, width, height);
		camera.aspect = ((float)width) / height;
	}

	void render_gui() {
		auto guiWindowWidth = 350;
		auto guiWindowHeight = 200;
		ImVec2 windowPosition = ImVec2(windowWidth - guiWindowWidth, 0);
		ImVec2 windowSize = ImVec2(guiWindowWidth, guiWindowHeight);

		ImGuiWindowFlags window_flags = 0;
		window_flags |= ImGuiWindowFlags_NoMove;
		window_flags |= ImGuiWindowFlags_NoResize;
		// window_flags |= ImGuiWindowFlags_NoCollapse;

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
		if (ImGui::SliderFloat("Noise Frequency", &perlinFreq, 0.001f, 0.1f)) generateTerrain();
		ImGui::Checkbox(" ", &brushFlag);
		ImGui::SameLine(0, 0);
		ImGui::SliderFloat("Brush Size", &brushSize, 1.0f, 10.0f);

		const char *items[] = {"NORMAL", "RANDOM"};
		static int selected_item = 0;
		static int previous_item = 0;

		if (ImGui::Combo("Biome Type", &selected_item, items, IM_ARRAYSIZE(items))) {
			if (selected_item == 0) biomeType = NORMAL;
			if (selected_item == 1) biomeType = RANDOM;
			if (previous_item != selected_item) {
				for (auto &vert : vertices)
					vert.color = heightmapColor(vert.pos, biomeType);

				gfx::upload(terrain, gfx::vertex_buffer, vertices);
			}
			previous_item = selected_item;
		}

		imguiWindowHovered = ImGui::IsWindowHovered() || ImGui::IsAnyItemHovered() || ImGui::IsMouseHoveringRect(windowPosition, ImVec2(windowPosition.x + windowSize.x, windowPosition.y + windowSize.y));

		ImGui::End();

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
	}

	void render_terrain() {
		gfx::layout<Vertex> basic_layout = {
			{0, &Vertex::pos},
			{1, &Vertex::color},
			{2, &Vertex::normal},
		};

		gfx::use(basic_shader);
		gfx::set_layout(terrain, basic_layout);
		pvm = camera.projection() * camera.view();
		nm = glm::mat4(1.0);
		eyePos = camera.position;
		gfx::draw(terrain);
	}

	void render_skybox() {
		gfx::layout<Vertex> skybox_layout = {
			{0, &Vertex::pos},
			{1, &Vertex::color},
		};
		gfx::use(skybox_shader);
		gfx::set_layout(skyBox, skybox_layout);
		gfx::draw(skyBox);
	}

	void render() {
		glClearColor(0, 0, 0, 1);
		glClearDepth(1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		render_skybox();
		render_terrain();
		render_gui();
	}

	void cameramove(float dt) {
		constexpr auto speed = 170;
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
		deltatime = dt;
		cameramove(dt);
		if (mouse.buttons[GLFW_MOUSE_BUTTON_LEFT] && !imguiWindowHovered && brushFlag) {
			terraformingOnClick(1);
		}
		if (mouse.buttons[GLFW_MOUSE_BUTTON_RIGHT] && !imguiWindowHovered && brushFlag) {
			terraformingOnClick(-1);
		}
		render();
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
	glfw::on(wnd, glfw::scroll, [&](glfw::window &window, double xoffset, double yoffset) { app.mousewheel(xoffset, yoffset); });

	app.init();

	glfw::run([&](float dt) {
		app.update(dt);

		glfw::swap_buffers(wnd);
		return glfw::is_closed(wnd);
	});

	glfw::imgui::uninitialize();

	return 0;
}
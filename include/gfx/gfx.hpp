#pragma once

#include <array>
#include <cmath>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>


#include <iostream>

#include <glad/glad.h>
#include <glm/ext.hpp>
#include <glm/glm.hpp>

namespace gfx {

	struct shader {
	  private:
		GLuint ID;

	  public:
		shader(GLuint ID = 0) noexcept : ID{ID} {}

		static shader create() { return shader(glCreateProgram()); }

		shader(const shader &) noexcept = delete;
		shader(shader &&o) noexcept : ID{std::exchange(o.ID, 0)} {}

		shader &operator=(const shader &) noexcept = delete;
		shader &operator=(shader &&o) noexcept {
			ID = std::exchange(o.ID, 0);
			return *this;
		}

		~shader() noexcept { glDeleteProgram(this->ID); }

		void attach(GLuint shaderType, const char *shaderCode) {
			GLuint shader = glCreateShader(shaderType);
			glShaderSource(shader, 1, std::addressof(shaderCode), nullptr);
			glCompileShader(shader);

			int success;
			glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
			if (!success) {
				char infoLog[512];
				glGetShaderInfoLog(shader, 512, nullptr, infoLog);
				std::cout << "ERROR::SHADER::COMPILATION_FAILED\n" << infoLog << std::endl;
				throw std::runtime_error("ERROR::SHADER::COMPILATION_FAILED\n");
			};

			glAttachShader(this->ID, shader);
			glDeleteShader(shader);
		}

		void link(std::string_view fragmentOutputName) {

			glBindFragDataLocation(this->ID, 0, std::data(fragmentOutputName));
			glLinkProgram(this->ID);

			int success;
			glGetProgramiv(this->ID, GL_LINK_STATUS, &success);

			if (!success) {
				char infoLog[512];
				glGetProgramInfoLog(this->ID, 512, nullptr, infoLog);
				std::cout << "ERROR::SHADER::COMPILATION_FAILED\n" << infoLog << std::endl;
				throw std::runtime_error("ERROR::SHADER::PROGRAM::LINKING_FAILED");
			}
		}

		void use() const noexcept { glUseProgram(this->ID); }

		constexpr explicit operator GLuint() const noexcept { return ID; }
	};

	template <typename T> struct uniform {};

	template <> struct uniform<glm::vec4> {

		GLuint location_ = 0;

		constexpr uniform() = default;
		uniform(const shader &shader_, const char *name_) noexcept {
			GLint res = glGetUniformLocation(static_cast<GLuint>(shader_), name_);
			if (res < 0)
				printf("uniform %s cannot be set\n", name_);
			else
				location_ = static_cast<GLuint>(res);
		}

		uniform &operator=(const uniform &value) = default;

		void operator=(glm::vec4 value) const noexcept { glUniform4fv(location_, 1, &value[0]); }
	};

	template <> struct uniform<glm::vec3> {

		GLuint location_ = 0;

		constexpr uniform() = default;
		uniform(const shader &shader_, const char *name_) noexcept {
			GLint res = glGetUniformLocation(static_cast<GLuint>(shader_), name_);
			if (res < 0)
				printf("uniform %s cannot be set\n", name_);
			else
				location_ = static_cast<GLuint>(res);
		}

		uniform &operator=(const uniform &value) = default;

		void operator=(glm::vec3 value) const noexcept { glUniform3fv(location_, 1, &value[0]); }
	};

	template <> struct uniform<glm::mat4> {

		GLuint location_ = 0;

		constexpr uniform() = default;
		uniform(const shader &shader_, const char *name_) noexcept {
			GLint res = glGetUniformLocation(static_cast<GLuint>(shader_), name_);
			if (res < 0)
				printf("uniform %s cannot be set\n", name_);
			else
				location_ = static_cast<GLuint>(res);
		}

		constexpr uniform(const uniform &value) : location_{value.location_} {}

		uniform &operator=(const uniform &value) noexcept {
			location_ = value.location_;
			return *this;
		}

		void operator=(const glm::mat4 &value) const noexcept { glUniformMatrix4fv(location_, 1, GL_FALSE, &value[0][0]); }
	};

	template <typename T> struct buffer {

		enum buffer_t { vertex, index };

		constexpr static auto max_buffer_count = 2;

		GLuint vertex_array_object_id = 0;
		GLsizei indexcount = 0;
		GLenum indextype = GL_UNSIGNED_SHORT;
		std::array<GLuint, max_buffer_count> vertex_buffer_object_ids = {0};
	};

	template <typename T> struct location {};

	template <typename T> struct is_glm_vector : std::false_type {};

	template <glm::length_t W, typename T, glm::qualifier Q> struct is_glm_vector<glm::vec<W, T, Q>> : std::true_type {};

	template <typename T> struct is_glm_matrix : std::false_type {};

	template <glm::length_t W, glm::length_t H, typename T, glm::qualifier Q> struct is_glm_matrix<glm::mat<W, H, T, Q>> : std::true_type {};

	template <typename T>
	concept glm_vector = is_glm_vector<T>::value;

	template <typename T>
	concept glm_matrix = is_glm_matrix<T>::value;

	template <typename T>
	concept glsl_basic_type = std::same_as<T, char> || std::same_as<T, unsigned char> || std::same_as<T, int> || std::same_as<T, unsigned int> || std::same_as<T, short> || std::same_as<T, unsigned short> || std::same_as<T, float>;

	template <typename T>
	concept glsl_type = glm_vector<T> || glm_matrix<T> || glsl_basic_type<T>;

	template <glsl_basic_type T> constexpr GLenum get_glsl_type();

	template <> constexpr GLenum get_glsl_type<float>() {
		return GL_FLOAT;
	}
	// template <> constexpr GLenum get_glsl_type<double>() { return GL_DOUBLE; }

	template <> constexpr GLenum get_glsl_type<char>() {
		return GL_BYTE;
	}
	template <> constexpr GLenum get_glsl_type<unsigned char>() {
		return GL_UNSIGNED_BYTE;
	}

	template <> constexpr GLenum get_glsl_type<short>() {
		return GL_SHORT;
	}
	template <> constexpr GLenum get_glsl_type<unsigned short>() {
		return GL_UNSIGNED_SHORT;
	}

	template <> constexpr GLenum get_glsl_type<int>() {
		return GL_INT;
	}
	template <> constexpr GLenum get_glsl_type<unsigned int>() {
		return GL_UNSIGNED_INT;
	}

	template <typename T> struct get_vec_length : std::integral_constant<glm::length_t, 0> {};

	template <glm::length_t W, typename T, glm::qualifier Q> struct get_vec_length<glm::vec<W, T, Q>> : std::integral_constant<glm::length_t, W> {};

	template <typename T> struct attribute_info;

	template <glsl_basic_type T> struct attribute_info<T> {
		constexpr static GLint count = 1;
		constexpr static GLint size = 1;
		constexpr static GLenum type = get_glsl_type<T::value_type>();
	};

	template <glm_vector T> struct attribute_info<T> {
		constexpr static GLint count = 1;
		constexpr static GLint size = get_vec_length<T>::value;
		constexpr static GLenum type = get_glsl_type<T::value_type>();
	};

	template <glm_matrix T> struct attribute_info<T> {
		constexpr static GLint count = get_vec_length<typename T::col_type>::value;
		constexpr static GLint size = get_vec_length<typename T::row_type>::value;
		constexpr static GLenum type = get_glsl_type<T::value_type>();
	};

	template <typename T> struct layout {

		struct attribute {

			template <glsl_type U> static auto member_to_offset(U T::*ptr) {
				T tmp;
				return static_cast<char *>(nullptr) + (reinterpret_cast<char *>(std::addressof(tmp.*ptr)) - reinterpret_cast<char *>(std::addressof(tmp)));
			}

			// clang-format off
			template <glsl_type U> 
			attribute(GLuint location, U T::*ptr) noexcept 
				: location(location)
				, size(attribute_info<U>::size)
				, type(attribute_info<U>::type)
				, count(attribute_info<U>::count)
				, offset(member_to_offset(ptr)) {}
			// clang-format on

			GLuint location;
			GLint size;
			GLenum type;
			GLint count;
			void *offset;
		};

		layout(const std::vector<attribute> &attributes) : attributes{attributes} {}
		layout(std::initializer_list<attribute> attributes) : attributes{attributes} {}

		GLsizei stride = sizeof(T);
		std::vector<attribute> attributes;
	};

	struct vertex_buffer_t {};
	struct index_buffer_t {};

	constexpr vertex_buffer_t vertex_buffer;
	constexpr index_buffer_t index_buffer;

	template <typename T, typename BuffType> void upload(gfx::buffer<BuffType> &buff, vertex_buffer_t, const std::vector<T> &data, GLenum usage = GL_STATIC_DRAW) noexcept {

		if (buff.vertex_array_object_id == 0) {
			glGenVertexArrays(1, &buff.vertex_array_object_id);
		}

		glBindVertexArray(buff.vertex_array_object_id);

		if (buff.vertex_buffer_object_ids[buffer<BuffType>::vertex] == 0) {
			glGenBuffers(1, &buff.vertex_buffer_object_ids[buffer<BuffType>::vertex]);
		}

		auto vbo = buff.vertex_buffer_object_ids[buffer<BuffType>::vertex];
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, std::size(data) * sizeof(T), std::data(data), usage);

		if (buff.vertex_buffer_object_ids[buffer<BuffType>::index] == 0) {
			buff.indexcount = std::size(data);
		}
	}

	template <typename T>
	concept opengl_index_type = std::same_as<T, unsigned char> || std::same_as<T, unsigned short> || std::same_as<T, unsigned int>;

	template <opengl_index_type T, typename BuffType> void upload(gfx::buffer<BuffType> &buff, index_buffer_t, const std::vector<T> &data, GLenum usage = GL_STATIC_DRAW) noexcept {

		if (buff.vertex_array_object_id == 0) {
			glGenVertexArrays(1, &buff.vertex_array_object_id);
		}

		glBindVertexArray(buff.vertex_array_object_id);

		if (buff.vertex_buffer_object_ids[buffer<BuffType>::index] == 0) {
			glGenBuffers(1, &buff.vertex_buffer_object_ids[buffer<BuffType>::index]);
		}

		auto ibo = buff.vertex_buffer_object_ids[buffer<BuffType>::index];
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, std::size(data) * sizeof(T), std::data(data), usage);

		buff.indexcount = std::size(data);
		buff.indextype = get_glsl_type<T>();
	}

	template <typename BuffType> void set_layout(gfx::buffer<BuffType> buff, gfx::layout<BuffType> layout) {

		if (buff.vertex_array_object_id == 0 || buff.vertex_buffer_object_ids[buffer<BuffType>::vertex] == 0) {
			throw std::runtime_error("Please upload vertex data in the buffer before use!");
		}

		glBindVertexArray(buff.vertex_array_object_id);
		for (auto &attr : layout.attributes) {
			for (GLint i = 0; i < attr.count; ++i) {
				glEnableVertexAttribArray(attr.location + i);
				glVertexAttribPointer(attr.location + i, attr.size, attr.type, GL_FALSE, layout.stride, attr.offset);
			}
		}
	}

	template <typename BuffType> void draw(gfx::buffer<BuffType> buff, GLenum mode = GL_TRIANGLES, GLsizei instancecount = 1) {

		if (buff.vertex_array_object_id == 0 || buff.vertex_buffer_object_ids[buffer<BuffType>::vertex] == 0) {
			return;
		}

		glBindVertexArray(buff.vertex_array_object_id);
		glDrawElementsInstanced(mode, buff.indexcount, buff.indextype, 0, instancecount);
	}

	template <typename T> void destroy(gfx::buffer<T> &buff) {
		if (buff.vertex_array_object_id != 0) {
			for (int i = 0; i < buff.max_buffer_count; ++i) {
				if (buff.vertex_buffer_object_ids[0] != 0) {
					glDeleteBuffers(1, &buff.vertex_buffer_object_ids[i]);
					buff.vertex_buffer_object_ids[i] = 0;
				}
			}
			glDeleteVertexArrays(1, &buff.vertex_array_object_id);
			buff.vertex_array_object_id = 0;
		}
	}

	glm::vec4 randVec4(float min = 0.0f, float max = 1.0f, float multiplier = 1.0f) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(min, max);
		return glm::vec4(dis(gen), dis(gen), dis(gen), dis(gen)) * multiplier;
	}

	void use(const gfx::shader &shader) {
		shader.use();
	}

	struct camera {

		constexpr static float YAW = -90.0f;
		constexpr static float PITCH = 0.0f;
		constexpr static float FOV = 45.0f;
		constexpr static float near = 0.1f;
		constexpr static float far = 200.0f;

		glm::vec3 position;
		glm::vec3 front;
		glm::vec3 up;
		glm::vec3 right;
		glm::vec3 worldup;

		float yaw;
		float pitch;
		float aspect = 800.0f / 600.0f;

		float fov;

		// clang-format off
		camera(glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f), float yaw = YAW, float pitch = PITCH) 
			: front(glm::vec3(0.0f, 0.0f, -1.0f))
			, fov(FOV) 
			, position(position)
			, worldup(up)
			, yaw(yaw)
			, pitch(pitch) {
			update_camera_vectors();
		}

		camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch) 
			: front(glm::vec3(0.0f, 0.0f, -1.0f))
			, fov(FOV)
			, position(posX, posY, posZ)
			, worldup(upX, upY, upZ)
			, yaw(yaw)
			, pitch(pitch) {
			update_camera_vectors();
		}
		// clang-format on

		void update_camera_vectors() {
			// calculate the new Front vector
			glm::vec3 front;
			front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
			front.y = sin(glm::radians(pitch));
			front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
			front = glm::normalize(front);
			right = glm::normalize(glm::cross(front, worldup));
			up = glm::normalize(glm::cross(right, front));
		}

		void lookat(const glm::vec3 &lookatpos) {
			front = glm::normalize(lookatpos - position);
			right = glm::normalize(glm::cross(front, worldup));
			up = glm::normalize(glm::cross(right, front));
			pitch = asin(front.y);
			yaw = glm::degrees(acos(front.x / cos(pitch)));
			pitch = glm::degrees(pitch);
		}

		glm::mat4 view() { return glm::lookAt(position, position + front, up); }
		glm::mat4 projection() { return glm::perspective(glm::radians(fov), aspect, near, far); }
	};
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

} // namespace gfx
#pragma once

#include <cassert>
#include <concepts>
#include <functional>
#include <stdexcept>

#include <GLFW/glfw3.h>

#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>
#include <imgui.h>

#include <glad/glad.h>

namespace glfw {

	struct context;

	context &instance();

	struct context {
	  private:
		context() {
			if (glfwInit() != GLFW_TRUE) {
				throw std::runtime_error("Couldn't initialize GLFW library!");
			}
		};

	  public:
		context(context &) = delete;
		context(context &&) = delete;

		context &operator=(context &) = delete;
		context &operator=(context &&) = delete;

		~context() noexcept { glfwTerminate(); }

		friend context &instance();
	};

	context &instance() {
		static context ctx{};
		return ctx;
	}

	// enum event_t { mousemove, resize, mousebutton, key, scroll, focus, close };

	struct event_mousemove_t {};
	struct event_resize_t {};
	struct event_mousebutton_t {};
	struct event_key_t {};
	struct event_scroll_t {};
	struct event_focus_t {};
	struct event_close_t {};

	constexpr event_mousemove_t mousemove{};
	constexpr event_resize_t resize{};
	constexpr event_mousebutton_t mousebutton{};
	constexpr event_key_t key{};
	constexpr event_scroll_t scroll{};
	constexpr event_focus_t focus{};
	constexpr event_close_t close{};

	struct window {

		window(int width = 640, int height = 480, std::string_view title = "GLFW Window") {
			window_ = glfwCreateWindow(width, height, std::data(title), nullptr, nullptr);
			glfwSetWindowUserPointer(window_, static_cast<void *>(this));
			init_callbacks();
		}

		window(window &) = delete;

		// clang-format off
		window(window &&other) noexcept
            : window_{ std::exchange(other.window_, nullptr) }
            , close{ std::move(other.close) }
            , focus{ std::move(other.focus) }
            , resize{ std::move(other.resize) }
            , mousebutton{ std::move(other.mousebutton) }
            , key{ std::move(other.key) }
            , mousemove{ std::move(other.mousemove) }
            , scroll{ std::move(other.scroll) }
		{}
        //clang-format on

		window &operator=(window &) = delete;
		window &operator=(window &&other) noexcept {
			window_ = std::exchange(other.window_, nullptr);
			close = std::move(other.close);
			focus = std::move(other.focus);
			resize = std::move(other.resize);
			mousebutton = std::move(other.mousebutton);
			key = std::move(other.key);
			mousemove = std::move(other.mousemove);
			scroll = std::move(other.scroll);
			return *this;
		}

		~window() noexcept {
			if (window_) glfwDestroyWindow(window_);
		}

		window &on(glfw::event_close_t, std::function<void(window &)> lambda) {
			close = std::move(lambda);
			return *this;
		}

		window &on(glfw::event_focus_t, std::function<void(window &, int)> lambda) {
			focus = std::move(lambda);
			return *this;
		}

		window &on(glfw::event_resize_t, std::function<void(window &, int, int)> lambda) {
			resize = std::move(lambda);
			return *this;
		}

		window &on(glfw::event_mousebutton_t, std::function<void(window &, int, int, int)> lambda) {
			mousebutton = std::move(lambda);
			return *this;
		}

		window &on(glfw::event_key_t, std::function<void(window &, int, int, int, int)> lambda) {
			key = std::move(lambda);
			return *this;
		}

		window &on(glfw::event_mousemove_t, std::function<void(window &, double, double)> lambda) {
			mousemove = std::move(lambda);
			return *this;
		}

		window &on(glfw::event_scroll_t, std::function<void(window &, double, double)> lambda) {
			scroll = std::move(lambda);
			return *this;
		}

        constexpr operator GLFWwindow*() const noexcept {
            return window_;
        }

	  private:
		GLFWwindow *window_ = nullptr;
		std::function<void(window &)> close;
		std::function<void(window &, int)> focus;
		std::function<void(window &, int, int)> resize;
		std::function<void(window &, int, int, int)> mousebutton;
		std::function<void(window &, int, int, int, int)> key;
		std::function<void(window &, double, double)> mousemove;
		std::function<void(window &, double, double)> scroll;

		void init_callbacks() {
			glfwSetFramebufferSizeCallback(window_, window::framebuffer_size_callback);
			glfwSetCursorPosCallback(window_, window::cursor_pos_callback);
			glfwSetMouseButtonCallback(window_, window::mouse_button_callback);
			glfwSetKeyCallback(window_, window::key_callback);
			glfwSetScrollCallback(window_, window::scroll_callback);
			glfwSetWindowFocusCallback(window_, window::focus_callback);
			glfwSetWindowCloseCallback(window_, window::close_callback);
		}

		static window &upcast(GLFWwindow *wnd) noexcept { return *static_cast<window *>(glfwGetWindowUserPointer(wnd)); }

		static void framebuffer_size_callback(GLFWwindow *wnd, int width, int height) {
			auto &self = upcast(wnd);
			if (self.resize) self.resize(self, width, height);
		}

		static void cursor_pos_callback(GLFWwindow *wnd, double xpos, double ypos) {
			auto &self = upcast(wnd);
			if (self.mousemove) self.mousemove(self, xpos, ypos);
		}

		static void mouse_button_callback(GLFWwindow *wnd, int button, int action, int mods) {
			auto &self = upcast(wnd);
			if (self.mousebutton) self.mousebutton(self, button, action, mods);
		}

		static void key_callback(GLFWwindow *wnd, int key, int scancode, int action, int mods) {
			auto &self = upcast(wnd);
			if (self.key) self.key(self, key, scancode, action, mods);
		}

		static void scroll_callback(GLFWwindow *wnd, double xoffset, double yoffset) {
			auto &self = upcast(wnd);
			if (self.scroll) self.scroll(self, xoffset, yoffset);
		}

		static void focus_callback(GLFWwindow *wnd, int focused) {
			auto &self = upcast(wnd);
			if (self.focus) self.focus(self, focused);
		}

		static void close_callback(GLFWwindow *wnd) {
			auto &self = upcast(wnd);
			if (self.close) self.close(self);
		}
	};

	struct mouse {
		float lastX;
		float lastY; 

    	float sensitivity = 0.1f;
		bool buttons[3] = {false};

	};

	struct keyboard {
		bool keys[349] = {false};
	};

    void on(glfw::window& wnd, glfw::event_close_t, std::function<void(window &)> lambda) {
        wnd.on(glfw::event_close_t{}, std::move(lambda));
    }

    void on(glfw::window& wnd, glfw::event_focus_t, std::function<void(window &, int)> lambda) {
        wnd.on(glfw::event_focus_t{}, std::move(lambda));
    }

    void on(glfw::window& wnd, glfw::event_resize_t, std::function<void(window &, int, int)> lambda) {
        wnd.on(glfw::event_resize_t{}, std::move(lambda));
    }

    void on(glfw::window& wnd, glfw::event_mousebutton_t, std::function<void(window &, int, int, int)> lambda) {
        wnd.on(glfw::event_mousebutton_t{}, std::move(lambda));
    }

    void on(glfw::window& wnd, glfw::event_key_t, std::function<void(window &, int, int, int, int)> lambda) {
        wnd.on(glfw::event_key_t{}, std::move(lambda));
    }

    void on(glfw::window& wnd, glfw::event_mousemove_t, std::function<void(window &, double, double)> lambda) {
        wnd.on(glfw::event_mousemove_t{}, std::move(lambda));
    }

    void on(glfw::window& wnd, glfw::event_scroll_t, std::function<void(window &, double, double)> lambda) {
        wnd.on(glfw::event_scroll_t{}, std::move(lambda));
    }

    void run(std::invocable<float> auto && lambda) {
        static float last = static_cast<float>(glfwGetTime());
        bool terminated = false;
        while (!terminated) {
            glfwPollEvents();
            float now = static_cast<float>(glfwGetTime());
            terminated = lambda(now - last);
            last = now;
        }
    }

    void swap_buffers(const glfw::window& wnd) {
        glfwSwapBuffers(wnd);
    }

    bool is_closed(const glfw::window& wnd) {
        return glfwWindowShouldClose(wnd);
    }

	namespace gl {

		void initialize() noexcept {
			glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
			glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
			glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
			glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, true);
		}
		void make_current(glfw::window& wnd) {
			glfwMakeContextCurrent(wnd);
		}

	} // namespace gl

	namespace imgui {

		void initialize(glfw::window& wnd) {
			IMGUI_CHECKVERSION();
			ImGui::CreateContext();
			ImGui_ImplGlfw_InitForOpenGL(wnd, true);
			ImGui_ImplOpenGL3_Init("#version 330");
		}

		void uninitialize() {
			ImGui_ImplOpenGL3_Shutdown();
			ImGui_ImplGlfw_Shutdown();
			ImGui::DestroyContext();
		}

	} 

	namespace glad {

		void initialize() {
			if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
				throw std::runtime_error("Failed to initialize GLAD");
			}
		}

	} 

} // namespace glfw
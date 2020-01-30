#ifndef VIEWER_H
#define VIEWER_H

#include "gl_includes.h"
#include "geometry.h"

#include <future>

struct viewer
{
	viewer(const int& nx, const int& ny, const int &ns, const int& num_channels);
	GLFWwindow* init();
<<<<<<< HEAD
	void background_thread(const std::shared_future<void> &future, GLFWwindow* window);
	void add_sample(const Vector2i &pixel, Vector3f sample);
	void save_and_destroy();

	// Stop rendering if to_exit is set to true by background_thread()
	bool to_exit;
=======
	static void background_thread(const std::shared_future<void> &future, GLubyte *out_image, GLFWwindow* window, int nx, int ny, bool &to_exit);
	void add_sample(const Vector2i &pixel, Vector3f sample);
	void save_and_destroy();

>>>>>>> 789dbd6... Refactor viewer into a separate class
	const int nx;
	const int ny;
	const int ns;
	const int num_channels;
	std::unique_ptr<GLubyte[]> out_image;
	std::unique_ptr<GLfloat[]> fout_image;
};

#endif
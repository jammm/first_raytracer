#ifndef VIEWER_H
#define VIEWER_H

#include "gl_includes.h"
#include "geometry.h"

#include <future>

struct viewer
{
	viewer(const int& nx, const int& ny, const int &ns, const int& num_channels);
	GLFWwindow* init();
	static void background_thread(const std::shared_future<void> &future, GLubyte *out_image, GLFWwindow* window, int nx, int ny, bool &to_exit);
	void add_sample(const Vector2i &pixel, Vector3f sample);
	void save_and_destroy();

	const int nx;
	const int ny;
	const int ns;
	const int num_channels;
	std::unique_ptr<GLubyte[]> out_image;
	std::unique_ptr<GLfloat[]> fout_image;
};

#endif
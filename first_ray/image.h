#ifndef IMAGE_H_
#define IMAGE_H_

#ifdef _WIN32
    #define GLEW_STATIC
    #define STBI_MSC_SECURE_CRT
#endif
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION

#include <string>
#include "stb_image.h"
#include "stb_image_write.h"


enum class stbi
{
    STBI_JPG,
    STBI_PNG,
    STBI_BMP
};

class image
{
    bool loaded_from_stbi;
public:
    image() {}
    image(const std::string &filename) : filename(filename), loaded_from_stbi(true)
    {
        data = stbi_load(filename.c_str(), &nx, &ny, &nn, 0);
    }

    image(unsigned char *data, const int &nx, const int &ny, const int &nn) : data(data), nx(nx), ny(ny), nn(nn),
																			  filename("out_test"), loaded_from_stbi(false)
	{}

    int save_image(stbi type)
    {
		std::string to_save;

        switch (type)
        {
            case stbi::STBI_JPG:
				
				to_save = (filename.find(".jpg") == std::string::npos) ? filename + ".jpg" : filename;
                stbi_write_bmp(to_save.c_str(), nx, ny, nn, (void *)data);
                break;
            case stbi::STBI_PNG:
				to_save = (filename.find(".png") == std::string::npos) ? filename + ".png" : filename;
                stbi_write_png(to_save.c_str(), nx, ny, nn, (void *)data, nx * nn);
                break;
            case stbi::STBI_BMP:
				to_save = (filename.find(".bmp") == std::string::npos) ? filename + ".bmp" : filename;
                stbi_write_jpg(to_save.c_str(), nx, ny, nn, (void *)data, 100);
                break;
            default:
                return -1;
        }

        return 0;
    }

    ~image()
    {
        if (loaded_from_stbi)
            stbi_image_free(data);
    }

    std::string filename;
    unsigned char *data;
    int nx;
    int ny;
    int nn;
};

#endif
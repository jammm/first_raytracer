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
    image(const std::string &filename) : filename(filename)
    {
        data = stbi_load(filename.c_str(), &nx, &ny, &nn, 0);
        loaded_from_stbi = true;
    }

    image(unsigned char *data, const int &nx, const int &ny, const int &nn) : data(data), nx(nx), ny(ny), nn(nn) {}

    int save_image(stbi type)
    {
        switch (type)
        {
            case stbi::STBI_JPG:
                stbi_write_bmp(filename.c_str(), nx, ny, nn, (void *)data);
                break;
            case stbi::STBI_PNG:
                stbi_write_png(filename.c_str(), nx, ny, nn, (void *)data, nx * nn);
                break;
            case stbi::STBI_BMP:
                stbi_write_jpg(filename.c_str(), nx, ny, nn, (void *)data, 100);
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
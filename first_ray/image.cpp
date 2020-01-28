#include "image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#ifdef _WIN32
    #define STBI_MSC_SECURE_CRT
#endif
#include "stb_image.h"
#include "stb_image_write.h"

image::image(const std::string &filename, const formats &format = formats::STBI_JPG) : loaded_from_stbi(true), filename(filename), type(format)
{
    switch (format)
    {
    case formats::STBI_HDR:
        dataf = stbi_loadf(filename.c_str(), &nx, &ny, &nn, 0);
        break;
    default:
        data = stbi_load(filename.c_str(), &nx, &ny, &nn, 0);
    }
}

int image::save_image(formats type)
{
    std::string to_save;

    // Flip output image vertically
    stbi_flip_vertically_on_write(true);

#ifdef __linux__
    filename.insert(0, "/");
    filename.insert(0, std::getenv("HOME"));
    std::cout<<"Filename is now "<<filename<<std::endl;
#endif

    switch (type)
    {
    case formats::STBI_JPG:

        to_save = (filename.find(".jpg") == std::string::npos) ? filename + ".jpg" : filename;
        stbi_write_bmp(to_save.c_str(), nx, ny, nn, (void *)data);
        break;
    case formats::STBI_PNG:
        to_save = (filename.find(".png") == std::string::npos) ? filename + ".png" : filename;
        stbi_write_png(to_save.c_str(), nx, ny, nn, (void *)data, nx * nn);
        break;
    case formats::STBI_BMP:
        to_save = (filename.find(".bmp") == std::string::npos) ? filename + ".bmp" : filename;
        stbi_write_jpg(to_save.c_str(), nx, ny, nn, (void *)data, 100);
        break;
    default:
        return -1;
    }

    return 0;
}

image::~image()
{
    if (loaded_from_stbi)
    {
        if (type == formats::STBI_HDR)
            stbi_image_free(dataf);
        else
            stbi_image_free(data);
    }
}
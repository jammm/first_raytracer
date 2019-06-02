#ifndef IMAGE_H_
#define IMAGE_H_

#include <vector>
#include <string>
#include <fstream>


enum class formats
{
    STBI_JPG,
    STBI_PNG,
    STBI_BMP,
};

class image
{
    bool loaded_from_stbi;
public:
    image() {}
    image(const std::string &filename);

    image(unsigned char *data, const int &nx, const int &ny, const int &nn) : data(data), nx(nx), ny(ny), nn(nn),
																			  filename("out_test"), loaded_from_stbi(false)
	{}

    int save_image(formats type);

    ~image();

    std::string filename;
    unsigned char *data;
    int nx;
    int ny;
    // number of channels
    int nn;
};

class image_pfm
{
public:
    image_pfm() {}
    image_pfm(float *data, const int &nx, const int &ny, const int &nn) : data(data, data + nx*ny*nn), nx(nx), ny(ny), nn(nn) {}
    image_pfm(const std::string &filename)
    {
        load_image(filename);
    }

    void load_image(std::string filename)
    {
        std::ifstream is(filename, std::ios::binary | std::ios::in);
        if (!is.is_open()) throw std::runtime_error("cannot find specified PFM file");

        // read pfm file header
        int &width = nx;
        int &height = ny;
        std::string pfm, minusOne;
        is >> pfm >> width >> height >> minusOne;
        char c = (char)is.get(); if (c == '\r') c = (char)is.get();

        // check pfm ill-conditions
        if (pfm != "PF" && pfm != "Pf") return;
        if (width < 0 || height < 0) return;

        nn = (pfm == "PF") ? 3 : 1;

        // allocate space for data
        data.reserve(width * height * nn);

        // read pfm file
        std::vector<float> row(nn * width, 0.0);
        for (int y = height; y > 0; y--)
        {
            is.read((char*)(&row[0]), sizeof(float) * nn * width);
            for (int x = 0; x < width; x++)
            {
                for (int c = 0; c < nn; c++)
                    data[(y*nx+x)*nn + c] = row[x*nn + c];
            }
        }

        // close
        is.close();
    }

    void save_image(std::string filename)
    {
        int &width = nx;
        int &height = ny;

        // write PFM header
        std::ofstream os(filename, std::ios::binary);
        if (!os.is_open()) throw std::runtime_error("cannot open file");
        os << "PF" << std::endl;
        os << width << " " << height << std::endl;
        os << "-1" << std::endl;

        // write data with row order flipped
        for (int y = 0; y < height; ++y)
            for (int x = 0; x < width; ++x)
                for (int c = 0; c < nn; ++c)
                {
                    const int index = (y*width + x)*nn + c;
                    const float vf = data[index];
                    os.write((char*)(&vf), sizeof(float));
                }

        os.close();
    }

    std::vector<float> data;
    int nx;
    int ny;
    // number of channels
    int nn;
};

#endif

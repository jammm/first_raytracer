#ifndef SAMPLER_H
#define SAMPLER_H

#include "geometry.h"
#include <thread>
#include <random>
#include <ctime>

struct sampler
{
    sampler(const int &seed) : engine(seed)
    {}

    sampler() : engine(std::time(0))
    {}

    inline double get1d()
    {
        return std::uniform_real_distribution<>(0.0, 1.0)(engine);
    }

    inline Vector2f get2d()
    {
        return Vector2f(std::uniform_real_distribution<>(0.0, 1.0)(engine), std::uniform_real_distribution<>(0.0, 1.0)(engine));
    }

    inline Vector3f get3d()
    {
        return Vector3f(std::uniform_real_distribution<>(0.0, 1.0)(engine),
                        std::uniform_real_distribution<>(0.0, 1.0)(engine),
                        std::uniform_real_distribution<>(0.0, 1.0)(engine));
    }

    std::mt19937 engine;
};

#endif
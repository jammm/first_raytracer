#ifndef UTIL_H
#define UTIL_H

#include <random>

// Generates the canonical uniform random variable ξ
inline float gen_cano_rand()
{
    thread_local static std::random_device seed_gen;
    thread_local static std::mt19937 engine(seed_gen());
    thread_local static std::uniform_real_distribution<> dist(0.0f, 1.0f);
    return dist(engine);
}

inline int box_x_compare(const void *a, const void *b)
{
    aabb box_left;
    aabb box_right;

    hitable *ah = *(hitable**)a;
    hitable *bh = *(hitable**)b;

    if (!ah->bounding_box(0.0f, 0.0f, box_left) || !bh->bounding_box(0.0f, 0.0f, box_right))
        std::cerr << "no bounding box in parallel_bvh_node constructor!\n";

    if (box_left.min.x() - box_right.min.x() < 0.0f)
        return -1;
    else
        return 1;
}

inline int box_y_compare(const void *a, const void *b)
{
    aabb box_left;
    aabb box_right;

    hitable *ah = *(hitable**)a;
    hitable *bh = *(hitable**)b;

    if (!ah->bounding_box(0.0f, 0.0f, box_left) || !bh->bounding_box(0.0f, 0.0f, box_right))
        std::cerr << "no bounding box in parallel_bvh_node constructor!\n";

    if (box_left.min.y() - box_right.min.y() < 0.0f)
        return -1;
    else
        return 1;
}

static inline float FromSrgb(const float v)
{
    if (v <= 0.04045) return v * (1.0 / 12.92);
    return std::pow((v + 0.055) * (1.0 / 1.055), 2.4);
}

inline int box_z_compare(const void *a, const void *b)
{
    aabb box_left;
    aabb box_right;

    hitable *ah = *(hitable**)a;
    hitable *bh = *(hitable**)b;

    if (!ah->bounding_box(0.0f, 0.0f, box_left) || !bh->bounding_box(0.0f, 0.0f, box_right))
        std::cerr << "no bounding box in parallel_bvh_node constructor!\n";

    if (box_left.min.z() - box_right.min.z() < 0.0f)
        return -1;
    else
        return 1;
}

#endif


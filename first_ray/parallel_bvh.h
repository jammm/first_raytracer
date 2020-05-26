#ifndef PARALLEL_BVH_H
#define PARALLEL_BVH_H

#include "hitable_list.h"
#include "util.h"
#include <taskflow/taskflow.hpp>

class parallel_bvh_node : public hitable_list
{
public:
    parallel_bvh_node() {}
    parallel_bvh_node(tf::SubflowBuilder &subflow, hitable **l, const int &n, const int &g_index, float time0, float time1);
    virtual bool hit(const ray &r, float tmin, float tmax, hit_record &rec) const;
    virtual bool bounding_box(float t0, float t1, aabb &b) const;

    static parallel_bvh_node *create_bvh(hitable **l, int n, float time0, float time1);

    // Members for each node
    std::unique_ptr<hitable> left;
    std::unique_ptr<hitable> right;
    aabb box;

    // Static members for the entire BVH
    static std::unique_ptr<aabb[]> g_boxes;
    static std::unique_ptr<float[]> g_left_area;
    static std::unique_ptr<float[]> g_right_area;
};

std::unique_ptr<aabb[]> parallel_bvh_node::g_boxes = nullptr;
std::unique_ptr<float[]> parallel_bvh_node::g_left_area = nullptr;
std::unique_ptr<float[]> parallel_bvh_node::g_right_area = nullptr;

bool parallel_bvh_node::bounding_box(float t0, float t1, aabb &b) const
{
    b = box;
    return true;
}

bool parallel_bvh_node::hit(const ray &r, float t_min, float t_max, hit_record &rec) const
{
    ++Scene::num_rays_traversed;
    // Martin Lamber's BVH improvement https://twitter.com/Peter_shirley/status/1105292977423900673
    if (box.hit(r, t_min, t_max))
    {
        ++Scene::num_ray_box_intersections;
        float ray_min_t = t_min;
        if (ray_min_t == EPSILON)
            ray_min_t *= std::max(std::max(std::max(std::abs(r.o[0]),
                std::abs(r.o[1])), std::abs(r.o[2])), EPSILON);

        if (ray_min_t > t_min) t_min = ray_min_t;

        if (left->hit(r, t_min, t_max, rec))
        {
            right->hit(r, t_min, rec.t, rec);
            return true;
        }
        else
        {
            return right->hit(r, t_min, t_max, rec);
        }
    }
    return false;
}

// Construct parallel_bvh using SAH method
parallel_bvh_node::parallel_bvh_node(tf::SubflowBuilder &subflow, hitable **l, const int &n, const int &g_index, float time0, float time1)
    : hitable_list(std::vector<hitable*>(l, l + n), n), left(nullptr), right(nullptr)
{
    ++Scene::num_bvh_nodes;
    if (!parallel_bvh_node::g_boxes)
    {
        parallel_bvh_node::g_boxes = std::make_unique<aabb[]>(n);
        parallel_bvh_node::g_left_area = std::make_unique<float[]>(n);
        parallel_bvh_node::g_right_area = std::make_unique<float[]>(n);
    }

    aabb* const boxes = g_boxes.get() + g_index;
    float* const left_area = g_left_area.get() + g_index;
    float* const right_area = g_right_area.get() + g_index;

    aabb main_box;
    l[0]->bounding_box(time0, time1, main_box);

    for (unsigned int i = 1; i < n; ++i)
    {
        aabb new_box;
        l[i]->bounding_box(time0, time1, new_box);
        main_box = surrounding_box(new_box, main_box);
    }

    int axis = main_box.longest_axis();
    if (axis == 0)
        qsort(l, n, sizeof(hitable *), box_x_compare);
    else if (axis == 1)
        qsort(l, n, sizeof(hitable *), box_y_compare);
    else
        qsort(l, n, sizeof(hitable *), box_z_compare);

    for (int i = 0; i < n; ++i)
        l[i]->bounding_box(time0, time1, boxes[i]);

    left_area[0] = boxes[0].area();
    aabb left_box = boxes[0];
    for (int i = 1; i < n - 1; ++i)
    {
        left_box = surrounding_box(left_box, boxes[i]);
        left_area[i] = left_box.area();
    }
    right_area[n - 1] = boxes[n - 1].area();
    aabb right_box = boxes[n - 1];
    for (int i = n - 2; i > 0; --i)
    {
        right_box = surrounding_box(right_box, boxes[i]);
        right_area[i] = right_box.area();
    }
    float min_SAH = FLT_MAX;
    int min_SAH_idx;
    for (int i = 0; i < n - 1; ++i)
    {
        float SAH = i * left_area[i] + (n - i - 1)*right_area[i + 1];
        if (SAH < min_SAH)
        {
            min_SAH_idx = i;
            min_SAH = SAH;
        }
    }

    if (min_SAH_idx == 0)
    {
        // reached a leaf node
        ++Scene::num_bvh_nodes;
        ++Scene::num_bvh_leaf_nodes;
        left.reset(l[0]);
    }
    else
    {
        subflow.emplace([&left = left , l, time0, time1, min_SAH_idx, g_index](tf::SubflowBuilder& subflow2)
        {
            left = std::make_unique<parallel_bvh_node>(subflow2, l, min_SAH_idx + 1, g_index, time0, time1);
        }).name("left");
    }
    if (min_SAH_idx == n - 2)
    {
        // reached another leaf node
        ++Scene::num_bvh_nodes;
        ++Scene::num_bvh_leaf_nodes;
        right.reset(l[min_SAH_idx + 1]);
    }
    else
    {
        subflow.emplace([&right = right, l, n, time0, time1, min_SAH_idx, g_index](tf::SubflowBuilder& subflow2) {
            right = std::make_unique<parallel_bvh_node>(subflow2, l + min_SAH_idx + 1, n - min_SAH_idx - 1, g_index + min_SAH_idx + 1, time0, time1);
        }).name("right");
    }

    subflow.detach();

    box = main_box;
}

parallel_bvh_node *parallel_bvh_node::create_bvh(hitable **l, int n, float time0, float time1)
{
    tf::Taskflow tf;
    parallel_bvh_node *root;
    auto task = tf.emplace([&](tf::SubflowBuilder &subflow)
    {
        root = new parallel_bvh_node(subflow, l, n, 0, time0, time1);
    }).name("root");

    tf.dispatch().get();
    //tf.dump_topologies(std::cout);
    
    return root;
}

#endif

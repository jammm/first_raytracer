#ifndef HITABLE_LIST_H_
#define HITABLE_LIST_H_

#include "hitable.h"
#include <assert.h>
#include <vector>

class hitable_list : public hitable
{
public:
    hitable_list() {}
    hitable_list(const std::vector<hitable *> &l , int n) { list = l; list_size = n; }

    bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override;
    bool bounding_box(double t0, double t1, aabb &b) const override;
    double pdf_direct_sampling(const hit_record &lrec, const Vector3f &to_light) const override;
    Vector3f sample_direct(hit_record &rec, const Vector3f &o, const Vector2f &sample) const override;
    inline hitable *operator[](const int &i) const { return list[i]; }

    /* Picks a random object from the list and returns its index (used for light sampling) */
    int pick_sample(const double &sample) const;
    int list_size;
    std::vector<hitable *> list;
};

#endif

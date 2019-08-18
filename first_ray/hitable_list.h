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

    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const;
    virtual bool bounding_box(float t0, float t1, aabb &b) const;
    virtual float pdf_direct_sampling(const hit_record &lrec, const Vector3f &to_light) const;
    virtual Vector3f sample_direct(hit_record &rec, const Vector3f &o) const;
    inline hitable *operator[](const int &i) const { return list[i]; }

    /* Picks a random object from the list and returns its index (used for light sampling) */
    int pick_sample() const;
    int list_size;
    std::vector<hitable *> list;
};

#endif

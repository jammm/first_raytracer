#ifndef BOX_H
#define BOX_H

#include "hitable_list.h"
#include "aarect.h"

class box : public hitable
{
public:
    box() {}
    box(const Vector3f &pmin, const Vector3f &pmax, material *mat);
    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const;
    virtual bool bounding_box(float t0, float t1, aabb &box) const
    {
        box = aabb(pmin, pmax);
        return true;
    }

    Vector3f pmin;
    Vector3f pmax;
    hitable *list_ptr;
};

box::box(const Vector3f &pmin, const Vector3f &pmax, material *mat) : pmin(pmin), pmax(pmax)
{
    hitable **list = new hitable*[6];
    list[0] = new xy_rect(pmin.x(), pmax.x(), pmin.y(), pmax.y(), pmax.z(), mat);
    list[1] = new flip_normals(new xy_rect(pmin.x(), pmax.x(), pmin.y(), pmax.y(), pmin.z(), mat));
    list[2] = new xz_rect(pmin.x(), pmax.x(), pmin.z(), pmax.z(), pmax.y(), mat);
    list[3] = new flip_normals(new xz_rect(pmin.x(), pmax.x(), pmin.z(), pmax.z(), pmin.y(), mat));
    list[4] = new yz_rect(pmin.y(), pmax.y(), pmin.z(), pmax.z(), pmax.x(), mat);
    list[5] = new flip_normals(new yz_rect(pmin.y(), pmax.y(), pmin.z(), pmax.z(), pmin.x(), mat));
    list_ptr = new hitable_list(list, 6);
}

bool box::hit(const ray &r, float t_min, float t_max, hit_record &rec) const
{
    return list_ptr->hit(r, t_min, t_max, rec);
}

#endif

#ifndef ONB_H_
#define ONB_H_

// Orthonormal basis generator
// Only works on directions
class onb
{
public:
    onb() {}
    inline Vector3f operator[](int i) const { return axis[i]; }
    Vector3f u() const { return axis[0]; }
    Vector3f v() const { return axis[1]; }
    Vector3f w() const { return axis[2]; }
    Vector3f local(const float &a, const float &b, const float &c) { return a*u() + b*v() + c*w(); }
    Vector3f local(const Vector3f &a) const { return a.x()*u() + a.y()*v() + a.z()*w(); }
    void onb::build_from_w(const Vector3f &n)
    {
        axis[2] = unit_vector(n);
        Vector3f a;
        if (fabs(w().x()) > 0.9)
            a = Vector3f(0, 1, 0);
        else
            a = Vector3f(1, 0, 0);

        axis[1] = unit_vector(cross(w(), a));
        axis[0] = cross(w(), v());
    }

    Vector3f axis[3];
};

#endif

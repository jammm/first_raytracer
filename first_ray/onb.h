#ifndef ONB_H_
#define ONB_H_

// Orthonormal basis generator
// Only works on directions
class onb
{
public:
    onb() {}
    onb(const Vector3f &w) { build_from_w(w); }
    inline Vector3f operator[](int i) const { return axis[i]; }
    Vector3f u() const { return axis[0]; }
    Vector3f v() const { return axis[1]; }
    Vector3f w() const { return axis[2]; }
    Vector3f local(const double &a, const double &b, const double &c) { return a*u() + b*v() + c*w(); }
    Vector3f fromLocal(const Vector3f &a) const { return a.x()*u() + a.y()*v() + a.z()*w(); }
    Vector3f toLocal(const Vector3f& a) const { return Vector3f(dot(a, axis[0]), dot(a, axis[1]), dot(a, axis[2])); }
    void build_from_w(const Vector3f &n)
    {
        axis[2] = n;
        if (std::abs(n[0]) > std::abs(n[1])) {
            double invLen = 1.0 / std::sqrt(n[0] * n[0] + n[2] * n[2]);
            axis[1] = Vector3f(n[2] * invLen, 0.0f, -n[0] * invLen);
        }
        else {
            double invLen = 1.0 / std::sqrt(n[1] * n[1] + n[2] * n[2]);
            axis[1] = Vector3f(0.0, n[2] * invLen, -n[1] * invLen);
        }
        axis[0] = cross(axis[1], axis[2]);
    }

    Vector3f axis[3];
};

#endif

#ifndef GEOMETRY_H
#define GEOMETRY_H


#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <stdlib.h>

template <typename T>
struct Vector3
{
    Vector3<T>() {}
    Vector3<T>(const float e1, const float e2, const float &e3) { e[0] = e1; e[1] = e2; e[2] = e3; }
    inline float x() const { return e[0]; }
    inline float y() const { return e[1]; }
    inline float z() const { return e[2]; }
    inline float r() const { return e[0]; }
    inline float g() const { return e[1]; }
    inline float b() const { return e[2]; }

    inline const Vector3<T>& operator+() const { return *this; }
    inline Vector3<T>operator-() const { return Vector3<T>(-e[0], -e[1], -e[2]); }
    inline float operator[](int i) const { return e[i]; }
    inline float& operator[](int i) { return e[i]; }

    inline Vector3<T>& operator+=(const Vector3<T> &v2);
    inline Vector3<T>& operator-=(const Vector3<T> &v2);
    inline Vector3<T>& operator*=(const Vector3<T> &v2);
    inline Vector3<T>& operator/=(const Vector3<T> &v2);
    inline Vector3<T>& operator*=(const T &t);
    inline Vector3<T>& operator/=(const T &t);

    inline T length() const
        { return static_cast<float>(sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2])); }
    inline T squared_length () const
        { return e[0] * e[0] + e[1] * e[1] + e[2] * e[2]; }
    inline void make_unit_vector();

    // e == 3 component vector in single precision floating point
    T e[3];
};

typedef Vector3<float> Vector3f;

template <typename T>
inline std::istream& operator >> (std::istream &is, Vector3<T> &t)
{
    is >> t.e[0] >> t.e[1] >> t.e[2];
    return is;
}

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const Vector3<T> &t)
{
    os << t.e[0] << " " << t.e[1] << " " << t.e[2];
    return os;
}

template <typename T>
inline void Vector3<T>::make_unit_vector()
{
    float k = 1.0f / sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
    e[0] *= k;
    e[1] *= k;
    e[2] *= k;
}

template <typename T>
inline Vector3<T> operator+(const Vector3<T>&v1, const Vector3<T>&v2)
{
    return Vector3<T>(v1.e[0] + v2.e[0], v1.e[1] + v2.e[1], v1.e[2] + v2.e[2]);
}

template <typename T>
inline Vector3<T> operator-(const Vector3<T>&v1, const Vector3<T>&v2)
{
    return Vector3<T>(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1], v1.e[2] - v2.e[2]);
}

template <typename T>
inline Vector3<T> operator*(const Vector3<T>&v1, const Vector3<T>&v2)
{
    return Vector3<T>(v1.e[0] * v2.e[0], v1.e[1] * v2.e[1], v1.e[2] * v2.e[2]);
}

template <typename T>
inline Vector3<T> operator/(const Vector3<T> &v1, const Vector3<T> &v2)
{
    return Vector3<T>(v1.e[0] / v2.e[0], v1.e[1] / v2.e[1], v1.e[2] / v2.e[2]);
}

template <typename T>
inline Vector3<T> operator*(const float &t, const Vector3<T>&v)
{
    return Vector3<T>(t*v.e[0], t*v.e[1], t*v.e[2]);
}

template <typename T>
inline Vector3<T>operator*(const Vector3<T>&v, const float &t)
{
    return t * v;
}

template <typename T>
inline Vector3<T>operator/(const Vector3<T>&v, float t)
{
    return Vector3<T>(v.e[0] / t, v.e[1] / t, v.e[2] / t);
}

template <typename T>
inline float dot(const Vector3<T>&v1, const Vector3<T>&v2)
{
    return v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2];
}

template <typename T>
inline Vector3<T>cross(const Vector3<T>&v1, const Vector3<T>&v2)
{
    return Vector3<T>(v1.e[1] * v2.e[2] - v1.e[2] * v2.e[1],
        -(v1.e[0] * v2.e[2] - v1.e[2] * v2.e[0]),
        v1.e[0] * v2.e[1] - v1.e[1] * v2.e[0]);
}

template <typename T>
inline Vector3<T>unit_vector(const Vector3<T>&v)
{
    return v / v.length();
}

template <typename T>
inline Vector3<T>& Vector3<T>::operator+=(const Vector3<T>& v2)
{
    e[0] += v2.e[0];
    e[1] += v2.e[1];
    e[2] += v2.e[2];
    return *this;
}

template <typename T>
inline Vector3<T>& Vector3<T>::operator-=(const Vector3<T>& v2)
{
    e[0] -= v2.e[0];
    e[1] -= v2.e[1];
    e[2] -= v2.e[2];
    return *this;
}

template <typename T>
inline Vector3<T>& Vector3<T>::operator*=(const Vector3<T>& v2)
{
    e[0] *= v2.e[0];
    e[1] *= v2.e[1];
    e[2] *= v2.e[2];
    return *this;
}

template <typename T>
inline Vector3<T>& Vector3<T>::operator/=(const Vector3<T>& v2)
{
    e[0] /= v2.e[0];
    e[1] /= v2.e[1];
    e[2] /= v2.e[2];
    return *this;
}

template <typename T>
inline Vector3<T>& Vector3<T>::operator*=(const T &t)
{
    e[0] *= t;
    e[1] *= t;
    e[2] *= t;
    return *this;
}

template <typename T>
inline Vector3<T>& Vector3<T>::operator/=(const T &t)
{
    float k = 1.0f / t;

    e[0] *= k;
    e[1] *= k;
    e[2] *= k;
    return *this;
}

#endif
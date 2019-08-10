#ifndef GEOMETRY_H
#define GEOMETRY_H


#define _USE_MATH_DEFINES
#include <cmath>
#if defined(_WIN32)
#include <math.h>
#endif
#include <iostream>
#include <stdlib.h>
#include <cassert>

template <typename T>

// Vector stuff
struct Vector3
{
    Vector3() {}
    Vector3(const float e1, const float e2, const float &e3) 
    { 
        e[0] = e1;
        e[1] = e2;
        e[2] = e3;

        assert(std::isfinite(e[0])
            && std::isfinite(e[1])
            && std::isfinite(e[2]));
    }

    template<typename U>
    explicit Vector3(const Vector3<U> &v)
    {
        e[0] = (T)v[0];
        e[1] = (T)v[1];
        e[2] = (T)v[2];

        assert(std::isfinite(e[0])
            && std::isfinite(e[1])
            && std::isfinite(e[2]));
    }

    inline float x() const { return e[0]; }
    inline float y() const { return e[1]; }
    inline float z() const { return e[2]; }
    inline float r() const { return e[0]; }
    inline float g() const { return e[1]; }
    inline float b() const { return e[2]; }

    inline const Vector3<T>& operator+() const { return *this; }
    inline Vector3 operator-() const { return Vector3<T>(-e[0], -e[1], -e[2]); }
    inline T operator[](int i) const { return e[i]; }
    inline T& operator[](int i) { return e[i]; }

    inline Vector3& operator+=(const Vector3<T> &v2);
    inline Vector3& operator-=(const Vector3<T> &v2);
    inline Vector3& operator*=(const Vector3<T> &v2);
    inline Vector3& operator/=(const Vector3<T> &v2);
    inline Vector3& operator*=(const T &t);
    inline Vector3& operator/=(const T &t);

    inline T length() const
        { return static_cast<float>(sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2])); }
    inline T squared_length () const
        { return e[0] * e[0] + e[1] * e[1] + e[2] * e[2]; }
    inline void make_unit_vector();

    // Vector3 data
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
    assert(std::isfinite(t));
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
    assert(std::isfinite(v2[0])
        && std::isfinite(v2[1])
        && std::isfinite(v2[2]));
    return v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2];
}

template <typename T>
inline Vector3<T>cross(const Vector3<T>&v1, const Vector3<T>&v2)
{
    return Vector3<T>(
        (v1.e[1] * v2.e[2]) - (v1.e[2] * v2.e[1]),
        (v1.e[2] * v2.e[0]) - (v1.e[0] * v2.e[2]),
        (v1.e[0] * v2.e[1]) - (v1.e[1] * v2.e[0])
    );
}

template <typename T>
inline Vector3<T>abs(const Vector3<T> &v)
{
    return Vector3<T>(abs(v.e[0]), abs(v.e[1]), abs(v.e[2]));
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

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2]));

    return *this;
}

template <typename T>
inline Vector3<T>& Vector3<T>::operator-=(const Vector3<T>& v2)
{
    e[0] -= v2.e[0];
    e[1] -= v2.e[1];
    e[2] -= v2.e[2];

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2]));

    return *this;
}

template <typename T>
inline Vector3<T>& Vector3<T>::operator*=(const Vector3<T>& v2)
{
    e[0] *= v2.e[0];
    e[1] *= v2.e[1];
    e[2] *= v2.e[2];

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2]));

    return *this;
}

template <typename T>
inline Vector3<T>& Vector3<T>::operator/=(const Vector3<T>& v2)
{
    e[0] /= v2.e[0];
    e[1] /= v2.e[1];
    e[2] /= v2.e[2];

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2]));

    return *this;
}

template <typename T>
inline Vector3<T>& Vector3<T>::operator*=(const T &t)
{
    e[0] *= t;
    e[1] *= t;
    e[2] *= t;

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2]));

    return *this;
}

template <typename T>
inline Vector3<T>& Vector3<T>::operator/=(const T &t)
{
    float k = 1.0f / t;

    e[0] *= k;
    e[1] *= k;
    e[2] *= k;

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2]));

    return *this;
}

// Point stuff
template <typename T>
struct Point2
{
    Point2() {}
    Point2(T _x, T _y) : x(_x), y(_y) {}

    template<typename U>
    explicit Point2(const Point2<U> &v)
    {
        x = (T)v.x;
        y = (T)v.y;
    }

    inline Point2<T>& operator*=(const Point2<T>& p2);
    inline T operator[](int i) const
    {
        if (i == 0) return x;
        return y;
    }
    inline T& operator[](int i)
    {
        if (i == 0) return x;
        return y;
    }


    // Point2 data
    T x, y;
};

template <typename T>
inline Point2<T> operator+(const Point2<T>&p1, const Point2<T>&p2)
{
    return Point2<T>(p1.x + p2.x, p1.y + p2.y);
}

template <typename T>
inline Point2<T> operator*(const Point2<T>&p1, const Point2<T>&p2)
{
    return Point2<T>(p1.x * p2.x, p1.y * p2.y);
}

template <typename T>
inline Point2<T>& Point2<T>::operator*=(const Point2<T>& p2)
{
    x *= p2.x;
    y *= p2.y;
    return *this;
}

template <typename T>
inline Point2<T> operator*(const float &t, const Point2<T>&p)
{
    return Point2<T>(t*p.x, t*p.y);
}

template <typename T>
inline Point2<T> operator*(const Point2<T>&p, const float &t)
{
    return t * p;
}

typedef Point2<float> Point2f;

#endif

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
#include <cstring>

namespace jam
{
    template <class T>
    constexpr bool isfinite(const T& val)
    {
        return val != std::numeric_limits<T>::quiet_NaN();
    }
}

template <typename T>
// Vector4 stuff
struct Vector4
{
    constexpr Vector4() = default;
    constexpr Vector4(const float e1, const float e2, const float e3, const float e4)
    {
        e[0] = e1;
        e[1] = e2;
        e[2] = e3;
        e[3] = e4;
    }

    template<typename U>
    explicit Vector4(const Vector4<U>& v)
    {
        e[0] = (T)v[0];
        e[1] = (T)v[1];
        e[2] = (T)v[2];
        e[3] = (T)v[3];

        assert(std::isfinite(e[0])
            && std::isfinite(e[1])
            && std::isfinite(e[2])
            && std::isfinite(e[3]));
    }

    inline float x() const { return e[0]; }
    inline float y() const { return e[1]; }
    inline float z() const { return e[2]; }
    inline float w() const { return e[3]; }
    inline float r() const { return e[0]; }
    inline float g() const { return e[1]; }
    inline float b() const { return e[2]; }
    inline float a() const { return e[3]; }

    inline const Vector4<T>& operator+() const { return *this; }
    inline Vector4 operator-() const { return Vector4<T>(-e[0], -e[1], -e[2], -e[3]); }
    inline T operator[](int i) const { return e[i]; }
    inline T& operator[](int i) { return e[i]; }

    inline Vector4& operator+=(const Vector4<T>& v2);
    inline Vector4& operator-=(const Vector4<T>& v2);
    inline Vector4& operator*=(const Vector4<T>& v2);
    inline Vector4& operator/=(const Vector4<T>& v2);
    inline Vector4& operator*=(const T& t);
    inline Vector4& operator/=(const T& t);

    inline T length() const
    {
        return static_cast<float>(sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2] + e[3] * e[3]));
    }
    inline T squared_length() const
    {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2] + e[3] * e[3];
    }
    inline void make_unit_vector();

    // Vector4 data
    T e[4] = { std::numeric_limits<T>::quiet_NaN(),
               std::numeric_limits<T>::quiet_NaN(),
               std::numeric_limits<T>::quiet_NaN(),
               std::numeric_limits<T>::quiet_NaN() };
};

typedef Vector4<float> Vector4f;

template <typename T>
inline std::istream& operator >> (std::istream& is, Vector4<T>& t)
{
    is >> t.e[0] >> t.e[1] >> t.e[2] >> t.e[3];
    return is;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Vector4<T>& t)
{
    os << t.e[0] << " " << t.e[1] << " " << t.e[2] << " " << t.e[3];
    return os;
}

template <typename T>
inline void Vector4<T>::make_unit_vector()
{
    T k = 1.0f / length();
    e[0] *= k;
    e[1] *= k;
    e[2] *= k;
    e[3] *= k;
}

template <typename T>
inline Vector4<T> operator+(const Vector4<T>& v1, const Vector4<T>& v2)
{
    return Vector4<T>(v1.e[0] + v2.e[0], v1.e[1] + v2.e[1], v1.e[2] + v2.e[2]);
}

template <typename T>
inline Vector4<T> operator-(const Vector4<T>& v1, const Vector4<T>& v2)
{
    return Vector4<T>(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1], v1.e[2] - v2.e[2], v1.e[3] - v2.e[3]);
}

template <typename T>
inline Vector4<T> operator*(const Vector4<T>& v1, const Vector4<T>& v2)
{
    return Vector4<T>(v1.e[0] * v2.e[0], v1.e[1] * v2.e[1], v1.e[2] * v2.e[2], v1.e[3] * v2.e[3]);
}

template <typename T>
inline Vector4<T> operator/(const Vector4<T>& v1, const Vector4<T>& v2)
{
    return Vector4<T>(v1.e[0] / v2.e[0], v1.e[1] / v2.e[1], v1.e[2] / v2.e[2], v1.e[3] / v2.e[3]);
}

template <typename T>
inline Vector4<T> operator*(const float& t, const Vector4<T>& v)
{
    return Vector4<T>(t * v.e[0], t * v.e[1], t * v.e[2], t * v.e[3]);
}

template <typename T>
inline Vector4<T>operator*(const Vector4<T>& v, const float& t)
{
    assert(std::isfinite(t));
    return t * v;
}

template <typename T>
inline Vector4<T>operator/(const Vector4<T>& v, float t)
{
    float denominator = 1.0f / t;
    return Vector4<T>(v.e[0] * t, v.e[1] * t, v.e[2] * t, v.e[3] * t);
}

template <typename T>
inline float dot(const Vector4<T>& v1, const Vector4<T>& v2)
{
    assert(std::isfinite(v2[0])
        && std::isfinite(v2[1])
        && std::isfinite(v2[2])
        && std::isfinite(v2[3]));
    return v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2] + v1.e[3] * v2.e[3];
}

template <typename T>
inline Vector4<T>abs(const Vector4<T>& v)
{
    return Vector4<T>(abs(v.e[0]), abs(v.e[1]), abs(v.e[2]), abs(v.e[3]));
}

template <typename T>
inline Vector4<T>unit_vector(const Vector4<T>& v)
{
    return v / v.length();
}

template <typename T>
inline Vector4<T>safe_sqrt(const Vector4<T>& v)
{
    Vector4f result;
    for (int i = 0; i < 4; ++i)
    {
        result[i] = std::max(static_cast<T>(0), std::sqrt(v[i]));
    }

    return result;
}

template <typename T>
inline Vector4<T>& Vector4<T>::operator+=(const Vector4<T>& v2)
{
    e[0] += v2.e[0];
    e[1] += v2.e[1];
    e[2] += v2.e[2];
    e[3] += v2.e[3];

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2])
        && std::isfinite(e[3]));

    return *this;
}

template <typename T>
inline Vector4<T>& Vector4<T>::operator-=(const Vector4<T>& v2)
{
    e[0] -= v2.e[0];
    e[1] -= v2.e[1];
    e[2] -= v2.e[2];
    e[3] -= v2.e[3];

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2])
        && std::isfinite(e[3]));

    return *this;
}

template <typename T>
inline Vector4<T>& Vector4<T>::operator*=(const Vector4<T>& v2)
{
    e[0] *= v2.e[0];
    e[1] *= v2.e[1];
    e[2] *= v2.e[2];
    e[3] *= v2.e[3];

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2])
        && std::isfinite(e[3]));

    return *this;
}

template <typename T>
inline Vector4<T>& Vector4<T>::operator/=(const Vector4<T>& v2)
{
    e[0] /= v2.e[0];
    e[1] /= v2.e[1];
    e[2] /= v2.e[2];
    e[3] /= v2.e[3];

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2])
        && std::isfinite(e[3]));

    return *this;
}

template <typename T>
inline Vector4<T>& Vector4<T>::operator*=(const T& t)
{
    e[0] *= t;
    e[1] *= t;
    e[2] *= t;
    e[3] *= t;

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2])
        && std::isfinite(e[3]));

    return *this;
}

template <typename T>
inline Vector4<T>& Vector4<T>::operator/=(const T& t)
{
    float k = 1.0f / t;

    e[0] *= k;
    e[1] *= k;
    e[2] *= k;
    e[3] *= k;

    assert(std::isfinite(e[0])
        && std::isfinite(e[1])
        && std::isfinite(e[2])
        && std::isfinite(e[3]));

    return *this;
}

template <typename T>
// Vector3 stuff
struct Vector3
{
    constexpr Vector3() = default;
    constexpr Vector3(const float e1, const float e2, const float e3) 
    { 
        e[0] = e1;
        e[1] = e2;
        e[2] = e3;
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
    T e[3] = {0, 0, 0};
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
inline Vector3<T>safe_sqrt(const Vector3<T>& v)
{
    Vector3f result;
    for (int i = 0; i < 3; ++i)
    {
        result[i] = std::max(static_cast<T>(0), std::sqrt(v[i]));
    }

    return result;
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

// 2D vector
template <typename T>
struct Vector2
{
    Vector2() {}
    Vector2(T _x, T _y) : x(_x), y(_y) {}

    template<typename U>
    explicit Vector2(const Vector2<U> &v)
    {
        x = (T)v.x;
        y = (T)v.y;
    }

    inline Vector2<T>& operator*=(const Vector2<T>& p2);
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


    // Vector2 data
    T x, y;
};

template <typename T>
inline Vector2<T> operator+(const Vector2<T>&p1, const Vector2<T>&p2)
{
    return Vector2<T>(p1.x + p2.x, p1.y + p2.y);
}

template <typename T>
inline Vector2<T> operator-(const Vector2<T>& p1, const Vector2<T>& p2)
{
    return Vector2<T>(p1.x - p2.x, p1.y - p2.y);
}

template <typename T>
inline Vector2<T> operator*(const Vector2<T>&p1, const Vector2<T>&p2)
{
    return Vector2<T>(p1.x * p2.x, p1.y * p2.y);
}

template <typename T>
inline Vector2<T>& Vector2<T>::operator*=(const Vector2<T>& p2)
{
    x *= p2.x;
    y *= p2.y;
    return *this;
}

template <typename T>
inline Vector2<T> operator*(const float &t, const Vector2<T>&p)
{
    return Vector2<T>(t*p[0], t*p[1]);
}

template <typename T>
inline Vector2<T> operator*(const Vector2<T>&p, const float &t)
{
    return t * p;
}

typedef Vector2<float> Vector2f;
typedef Vector2<int> Vector2i;

// Matrix class taken from Mitsuba
template <int M, int N, typename T> struct Matrix {
    T m[M][N];

    Matrix() 
    {
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                m[i][j] = std::numeric_limits<double>::quiet_NaN();
    }

    explicit inline Matrix(T value) 
    {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] = value;
    }

    explicit inline Matrix(const T _m[M][N])
    {
        std::memcpy(m, _m, sizeof(T) * M * N);
    }

    explicit inline Matrix(const T _m[M * N])
    {
        std::memcpy(m, _m, sizeof(T) * M * N);
    }

    inline Matrix(const Matrix& mtx) 
    {
        std::memcpy(m, mtx.m, sizeof(T) * M * N);
    }

    void set_identity()
    {
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                m[i][j] = (i == j) ? 1.0f : 0.0f;
    }

    inline void transpose(Matrix<N, M, T> &target) const {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                target.m[i][j] = m[j][i];
    }

    void set_zero()
    {
        memset(m, 0, sizeof(T) * M * N);
    }

    bool invert(Matrix& target) const;
    inline T &operator()(int i, int j) { return m[i][j]; }
    inline const T & operator()(int i, int j) const { return m[i][j]; }

    inline bool operator==(const Matrix &mat) const 
    {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                if (m[i][j] != mat.m[i][j])
                    return false;
        return true;
    }

    inline bool operator!=(const Matrix &mat) const 
    {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                if (m[i][j] != mat.m[i][j])
                    return true;
        return false;
    }

    inline Matrix &operator=(const Matrix &mat) 
    {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] = mat.m[i][j];
        return *this;
    }

    inline Matrix operator+(const Matrix &mat) const 
    {
        Matrix result;
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                result.m[i][j] = m[i][j] + mat.m[i][j];
        return result;
    }

    inline Matrix operator+(T value) const 
    {
        Matrix result;
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                result.m[i][j] = m[i][j] + value;
        return result;
    }

    inline const Matrix &operator+=(const Matrix &mat) 
    {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] += mat.m[i][j];
        return *this;
    }

    inline const Matrix &operator+=(T value) 
    {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] += value;
        return *this;
    }

    inline Matrix operator-(const Matrix &mat) const 
    {
        Matrix result;
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                result.m[i][j] = m[i][j] - mat.m[i][j];
        return result;
    }

    inline Matrix operator-(T value) const 
    {
        Matrix result;
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                result.m[i][j] = m[i][j] - value;
        return result;
    }

    inline const Matrix &operator-=(const Matrix &mat) 
    {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] -= mat.m[i][j];
        return *this;
    }

    inline const Matrix &operator-(T value) 
    {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] -= value;
        return *this;
    }

    inline const Matrix &operator-=(T value) 
    {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] -= value;
        return *this;
    }

    inline Matrix operator-() const 
    {
        Matrix result;
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                result.m[i][j] = -m[i][j];
        return result;
    }

    inline Matrix operator*(T value) const 
    {
        Matrix result;
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                result.m[i][j] = m[i][j]*value;
        return result;
    }

    inline const Matrix& operator*=(T value) 
    {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] *= value;
        return *this;
    }

    inline Matrix operator/(T value) const 
    {
        Matrix result;
        float denominator = 1/value;
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                result.m[i][j] = m[i][j]*denominator;
        return result;
    }

    inline const Matrix& operator/=(T value) 
    {
        float denominator = 1/value;
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] *= denominator;
        return *this;
    }

    inline const Matrix &operator*=(const Matrix &mat) 
    {
        static_assert(M == N);
        Matrix temp = *this * mat;
        *this = temp;
        return *this;
    }

    inline Vector3f normal_transform(const Vector3f &v) const
    {

        float x = m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2];
        float y = m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2];
        float z = m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2];

        return Vector3f(x, y, z);
    }
};

struct Matrix4x4 : public Matrix<4, 4, float> 
{
    inline Matrix4x4() { }

    explicit inline Matrix4x4(float value) : Matrix<4, 4, float>(value) { }

    explicit inline Matrix4x4(const float _m[4][4]) : Matrix<4, 4, float>(_m) { }

    explicit inline Matrix4x4(const float _m[16]) : Matrix<4, 4, float>(_m) { }

    explicit inline Matrix4x4(const Vector4f &v1, const Vector4f &v2, const Vector4f &v3, const Vector4f &v4) 
    {
        m[0][0] = v1[0]; m[0][1] = v2[0]; m[0][2] = v3[0]; m[0][3] = v4[0];
        m[1][0] = v1[1]; m[1][1] = v2[1]; m[1][2] = v3[1]; m[1][3] = v4[1];
        m[2][0] = v1[2]; m[2][1] = v2[2]; m[2][2] = v3[2]; m[2][3] = v4[2];
        m[3][0] = v1[3]; m[3][1] = v2[3]; m[3][2] = v3[3]; m[3][3] = v4[3];
    }

    inline Matrix4x4(const Matrix<4, 4, float> &mtx) : Matrix<4, 4, float>(mtx) {}

    inline Matrix4x4(
        float a00, float a01, float a02, float a03,
        float a10, float a11, float a12, float a13,
        float a20, float a21, float a22, float a23,
        float a30, float a31, float a32, float a33) 
    {
        m[0][0] = a00; m[0][1] = a01; m[0][2] = a02; m[0][3] = a03;
        m[1][0] = a10; m[1][1] = a11; m[1][2] = a12; m[1][3] = a13;
        m[2][0] = a20; m[2][1] = a21; m[2][2] = a22; m[2][3] = a23;
        m[3][0] = a30; m[3][1] = a31; m[3][2] = a32; m[3][3] = a33;
    }

    inline float det3x3() const
    {
        return ((m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]))
              - (m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]))
              + (m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])));
    }

    inline Vector4f operator*(const Vector4f &v) const
    {
        return Vector4f(
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2] + m[0][3] * v[3],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2] + m[1][3] * v[3],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2] + m[2][3] * v[3],
            m[3][0] * v[0] + m[3][1] * v[1] + m[3][2] * v[2] + m[3][3] * v[3]
        );
    }

    inline Vector3f operator*(const Vector3f &v) const
    {

        float x = m[0][0] * v[0] + m[0][1] * v[1]
            + m[0][2] * v[2] + m[0][3];
        float y = m[1][0] * v[0] + m[1][1] * v[1]
            + m[1][2] * v[2] + m[1][3];
        float z = m[2][0] * v[0] + m[2][1] * v[1]
            + m[2][2] * v[2] + m[2][3];
        float w = m[3][0] * v[0] + m[3][1] * v[1]
            + m[3][2] * v[2] + m[3][3];

        if (w == 1.0f)
            return Vector3f(x, y, z);
        
        return Vector3f(x, y, z) / w;
    }

    inline Matrix4x4 operator*(float value) const 
    {
        Matrix4x4 result;
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                result.m[i][j] = m[i][j]*value;
        return result;
    }

    inline Matrix4x4 &operator=(const Matrix<4, 4, float> &mat) 
    {
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                m[i][j] = mat.m[i][j];
        return *this;
    }

    inline Vector4f row(int i) const 
    {
        return Vector4f(
            m[i][0], m[i][1], m[i][2], m[i][3]
        );
    }

    inline Vector4f col(int i) const 
    {
        return Vector4f(
            m[0][i], m[1][i], m[2][i], m[3][i]
        );
    }
};

template <typename T, int M1, int N1, int M2, int N2> inline Matrix<M1, N2, T>
        operator*(const Matrix<M1, N1, T> &mat1, const Matrix<M2, N2, T> &mat2) {
    static_assert(N1 == M2);
    Matrix<M1, N2, T> result;
    for (int i=0; i<M1; ++i) {
        for (int j=0; j<N2; ++j) {
            T sum = 0;
            for (int k=0; k<N1; ++k)
                sum += mat1.m[i][k] * mat2.m[k][j];
            result.m[i][j] = sum;
        }
    }
    return result;
}

// Matrix inversion taken from mitsuba
template <int M, int N, typename T> bool Matrix<M, N, T>::invert(Matrix &target) const {
    static_assert(M == N);

    int indxc[N], indxr[N];
    int ipiv[N];
    memset(ipiv, 0, sizeof(int)*N);
    std::memcpy(target.m, m, M*N*sizeof(T));

    for (int i = 0; i < N; i++) {
        int irow = -1, icol = -1;
        T big = 0;
        for (int j = 0; j < N; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < N; k++) {
                    if (ipiv[k] == 0) {
                        if (std::abs(target.m[j][k]) >= big) {
                            big = std::abs(target.m[j][k]);
                            irow = j;
                            icol = k;
                        }
                    } else if (ipiv[k] > 1) {
                        return false;
                    }
                }
            }
        }
        ++ipiv[icol];
        if (irow != icol) {
            for (int k = 0; k < N; ++k)
                std::swap(target.m[irow][k], target.m[icol][k]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (target.m[icol][icol] == 0)
            return false;
        T pivinv = 1.f / target.m[icol][icol];
        target.m[icol][icol] = 1.f;
        for (int j = 0; j < N; j++)
            target.m[icol][j] *= pivinv;
        for (int j = 0; j < N; j++) {
            if (j != icol) {
                T save = target.m[j][icol];
                target.m[j][icol] = 0;
                for (int k = 0; k < N; k++)
                    target.m[j][k] -= target.m[icol][k]*save;
            }
        }
    }
    for (int j = N-1; j >= 0; j--) {
        if (indxr[j] != indxc[j]) {
            for (int k = 0; k < N; k++)
                std::swap(target.m[k][indxr[j]], target.m[k][indxc[j]]);
        }
    }
    return true;
}

template <typename T, int M, int N> inline Matrix<M, N, T> operator*(T f, const Matrix<M, N, T> &m) {
    return m*f;
}

#endif

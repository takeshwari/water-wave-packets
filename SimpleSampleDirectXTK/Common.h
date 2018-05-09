#pragma once

#include <Eigen/Core>

#include <tinyformat.h>

#include <tbb/parallel_for.h>

#include <cmath>
#include <exception>
#include <string>
#include <vector>

// Enable TBB parallelization
#define USE_TBB 1

namespace pbs {

// Types ----------------------------------------------------------------------

template <typename Scalar, int Dimension>  struct TVector;
template <typename Vector>                 struct TBox;

typedef TVector<float, 1>       Vector1f;
typedef TVector<float, 2>       Vector2f;
typedef TVector<float, 3>       Vector3f;
typedef TVector<float, 4>       Vector4f;
typedef TVector<double, 1>      Vector1d;
typedef TVector<double, 2>      Vector2d;
typedef TVector<double, 3>      Vector3d;
typedef TVector<double, 4>      Vector4d;
typedef TVector<int, 1>         Vector1i;
typedef TVector<int, 2>         Vector2i;
typedef TVector<int, 3>         Vector3i;
typedef TVector<int, 4>         Vector4i;
typedef TVector<uint32_t, 1>    Vector1u;
typedef TVector<uint32_t, 2>    Vector2u;
typedef TVector<uint32_t, 3>    Vector3u;
typedef TVector<uint32_t, 4>    Vector4u;
typedef TBox<Vector1f>          Box1f;
typedef TBox<Vector2f>          Box2f;
typedef TBox<Vector3f>          Box3f;
typedef TBox<Vector4f>          Box4f;
typedef TBox<Vector1d>          Box1d;
typedef TBox<Vector2d>          Box2d;
typedef TBox<Vector3d>          Box3d;
typedef TBox<Vector4d>          Box4d;
typedef TBox<Vector1i>          Box1i;
typedef TBox<Vector2i>          Box2i;
typedef TBox<Vector3i>          Box3i;
typedef TBox<Vector4i>          Box4i;

typedef Eigen::Matrix<float,    Eigen::Dynamic, Eigen::Dynamic> MatrixXf;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXu;

// Math utilities -------------------------------------------------------------

template<typename T>
static constexpr T sqr(T x) { return x*x; }

template<typename T>
static constexpr T cube(T x) { return x*x*x; }
/*
template<typename T>
static inline T clamp(T x, T lo, T hi)  {
    return std::max(lo, std::min(hi, x));
}
*/
template<typename S, typename T>
static inline T lerp(S t, const T &a, const T &b) {
    return (S(1) - t) * a + t * b;
}
/*
static inline float unitToRange(float x, float lo, float hi) {
    return lo + clamp(x, 0.f, 1.f) * (hi - lo);
}

static inline float rangeToUnit(float x, float lo, float hi) {
    return clamp((x - lo) / (hi - lo), 0.f, 1.f);
}
*/
static inline uint32_t nextPowerOfTwo(uint32_t x) {
    x--;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x++;
    return x;
}

// String utilities -----------------------------------------------------------

std::vector<std::string> tokenize(const std::string &s, const std::string &delim = ", ", bool includeEmpty = false);

std::string toLower(const std::string &value);
bool toBool(const std::string &str);
int toInt(const std::string &str);
unsigned int toUInt(const std::string &str);
float toFloat(const std::string &str);

// Convert a time value in milliseconds into a human-readable string
std::string timeString(double time, bool precise = false);

// Convert a memory amount in bytes into a human-readable string
std::string memString(size_t size, bool precise = false);

/// Indent a string by the specified number of spaces
std::string indent(const std::string &string, int amount = 2);

// Threading

// iterate i=0..count-1 calling func(i)
template<typename Func>
inline void parallelFor(size_t count, Func func) {
#if USE_TBB
    tbb::parallel_for(0ul, count, 1ul, [func] (size_t i) { func(i); });
#else
    for (size_t i = 0; i < count; ++i) { func(i); }
#endif
}

// Debugging ------------------------------------------------------------------

class Exception : public std::runtime_error {
public:
    template<typename... Args>
    Exception(const char *fmt, const Args &... args) :
        std::runtime_error(tfm::format(fmt, args...))
    {}
};

template<typename... Args>
static inline void DBG(const char *fmt, const Args &... args) {
    std::cout << tfm::format(fmt, args...) << std::endl;
}

template<typename... Args>
static inline void handleAssert(bool cond, const char *fmt, const Args &... args) {
    if (!cond) {
        std::string msg = tfm::format(fmt, args...);
        std::cout << msg << std::endl;
        throw std::runtime_error(msg);
    }
}

#define ASSERT(_cond_, _fmt_, _args_) \
handleAssert(_cond_, "%s(%d) ASSERT: " _fmt_, __FILE__, __LINE__, ##_args_)


} // namespace pbs

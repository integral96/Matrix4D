#pragma once

#define TBB_PREVIEW_BLOCKED_RANGE_ND 1

#include <type_traits>

#include <boost/multi_array.hpp>
#include <boost/mpl/int.hpp>
#include <boost/proto/core.hpp>
#include <boost/proto/context.hpp>
#include <boost/random.hpp>
#include <boost/noncopyable.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_rangeNd.h>

namespace mpl = boost::mpl;
namespace proto = boost::proto;

namespace _spatial {

/////////MATRIX1D
template<typename T>
struct matrix1_ : boost::multi_array<T, 1> {
    typedef boost::multi_array<T, 1> array_type;
    typedef boost::multi_array_types::index_range range;
    typedef typename array_type::template array_view<1>::type sub_matrix_view1D;
    typedef T value_type;
public:
    constexpr matrix1_(const std::array<size_t, 1>& shape) : array_type(shape) {}
    constexpr matrix1_(size_t N) : array_type({boost::extents[N]}) {}

    constexpr void is_Matrix() {}

    inline T& operator () (size_t i) {
        return (*this)[i];
    }
    inline T const& operator () (size_t i) const {
        return (*this)[i];
    }
};

/////////MATRIX2D
template<typename T>
struct matrix2_ : boost::multi_array<T, 2> {
    typedef boost::multi_array<T, 2> array_type;
    typedef boost::multi_array_types::index_range range;
    typedef typename array_type::template array_view<1>::type sub_matrix_view1D;
    typedef typename array_type::template array_view<2>::type sub_matrix_view2D;
    typedef T value_type;
public:
    constexpr matrix2_(const std::array<size_t, 2>& shape) : array_type(shape) {}
    constexpr matrix2_(size_t N, size_t M) : array_type({boost::extents[N][M]}) {}

    constexpr void is_Matrix() {}

    inline T& operator () (size_t i, size_t j) {
        return (*this)[i][j];
    }
    inline T const& operator () (size_t i, size_t j) const {
        return (*this)[i][j];
    }
};

/////////MATRIX3D
template<typename T>
struct matrix3_ : boost::multi_array<T, 3> {
    typedef boost::multi_array<T, 3> array_type;
    typedef boost::multi_array_types::index_range range;
    typedef typename array_type::template array_view<1>::type sub_matrix_view1D;
    typedef typename array_type::template array_view<2>::type sub_matrix_view2D;
    typedef typename array_type::template array_view<3>::type sub_matrix_view3D;
    typedef T value_type;
public:
    constexpr matrix3_(const std::array<size_t, 3>& shape) : array_type(shape) {}
    constexpr matrix3_(size_t N, size_t M, size_t K) : array_type({boost::extents[N][M][K]}) {}

    constexpr void is_Matrix() {}

    inline T& operator () (size_t i, size_t j, size_t k) {
        return (*this)[i][j][k];
    }
    inline T const& operator () (size_t i, size_t j, size_t k) const {
        return (*this)[i][j][k];
    }
};

/////////MATRIX4D
template<typename T>
struct matrix4_ : boost::multi_array<T, 4> {
    typedef boost::multi_array<T, 4> array_type;
    typedef boost::multi_array_types::index_range range;
    typedef typename array_type::template array_view<1>::type sub_matrix_view1D;
    typedef typename array_type::template array_view<2>::type sub_matrix_view2D;
    typedef typename array_type::template array_view<3>::type sub_matrix_view3D;
    typedef typename array_type::template array_view<4>::type sub_matrix_view4D;
    typedef T value_type;
    array_type MTRX;
public:
    constexpr matrix4_(const std::array<size_t, 4>& shape) : array_type(shape) {}
    constexpr matrix4_(size_t N, size_t M, size_t K, size_t L) : array_type({boost::extents[N][M][K][L]}) {}

    constexpr void is_Matrix() {}

    inline T& operator () (size_t i, size_t j, size_t k, size_t l) {
        return (*this)[i][j][k][l];
    }
    inline T const& operator () (size_t i, size_t j, size_t k, size_t l) const {
        return (*this)[i][j][k][l];
    }
};

template<typename T>
struct Matrix1D;
template<typename T>
struct Matrix2D;
template<typename T>
struct Matrix3D;
template<typename T>
struct Matrix4D;

}

template<typename Matrix, class = void>
struct IsMatrix : mpl::false_ {};
template<typename Matrix>
struct IsMatrix<Matrix, std::void_t<decltype (std::declval<Matrix&>().is_Matrix())> > : mpl::true_  {};

template<typename T, class Matrix>
struct IsMatrix1D : mpl::false_ {};
template<typename T>
struct IsMatrix1D<T, _spatial::Matrix1D<T>> : mpl::true_ {};

template<typename T, class Matrix>
struct IsMatrix2D : mpl::false_ {};
template<typename T>
struct IsMatrix2D<T, _spatial::Matrix2D<T>> : mpl::true_ {};

template<typename T, class Matrix>
struct IsMatrix3D : mpl::false_ {};
template<typename T>
struct IsMatrix3D<T, _spatial::Matrix3D<T>> : mpl::true_ {};

template<typename T, class Matrix>
struct IsMatrix4D : mpl::false_ {};
template<typename T>
struct IsMatrix4D<T, _spatial::Matrix4D<T>> : mpl::true_ {};

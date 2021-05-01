#pragma once

#include <include/Matrix1D.hpp>
#include <include/Matrix2D.hpp>
#include <include/Matrix3D.hpp>
#include <include/Matrix4D.hpp>

template<typename T>
using sec_matrix = mpl::vector<_spatial::Matrix1D<T>, _spatial::Matrix2D<T>, _spatial::Matrix3D<T>, _spatial::Matrix4D<T>>;

template<size_t N, typename T, typename = std::enable_if_t<(N > 0)>>
using MATRIX = typename  mpl::at_c<sec_matrix<T>, N - 1>::type;

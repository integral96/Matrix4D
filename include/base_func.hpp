#pragma once

#include <iostream>

#include <boost/mpl/bool.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/type_traits/enable_if.hpp>
#include <boost/multi_array.hpp>
#include <boost/core/demangle.hpp>
#include <boost/core/typeinfo.hpp>


/*!
 * struct meta_loop
 */
template <size_t N, size_t I, class Closure>
typename boost::enable_if_t<(I == N)> is_meta_loop(Closure& closure) {}

template <size_t N, size_t I, class Closure>
typename boost::enable_if_t<(I < N)> is_meta_loop(Closure& closure) {
    closure.template apply<I>();
    is_meta_loop<N, I + 1>(closure);
}
template <size_t N, class Closure>
void meta_loop(Closure& closure) {
    is_meta_loop<N, 0>(closure);
}
template <size_t N, class Closure>
void meta_loopUV(Closure& closure) {
    is_meta_loop<N, 1>(closure);
}

/*!
 * struct abstract_sum
 */
template<class Closure>
struct abstract_sum_closures {
    typedef typename Closure::value_type value_type;
    abstract_sum_closures(Closure &closure) :  closure(closure), result(value_type()){}

    template<unsigned I>
    void apply(){
        result += closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};

template<unsigned N, class Closure>
typename Closure::value_type abstract_sums(Closure &closure) {
    abstract_sum_closures<Closure> my_closure(closure);
    meta_loop<N>(my_closure);
    return my_closure.result;
}

/*!
 * struct abstract_subtract
 */
template<class Closure>
struct abstract_subtract_closures {
    typedef typename Closure::value_type value_type;
    abstract_subtract_closures(Closure &closure) :  closure(closure), result(value_type()){}

    template<unsigned I>
    void apply(){
        result -= closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};

template<unsigned N, class Closure>
typename Closure::value_type abstract_subtract(Closure &closure) {
    abstract_subtract_closures<Closure> my_closure(closure);
    meta_loop<N>(my_closure);
    return my_closure.result;
}

/*!
 * struct abstract_divide
 */
template<class Closure>
struct abstract_divide_closures {
    typedef typename Closure::value_type value_type;
    abstract_divide_closures(Closure &closure) :  closure(closure), result(value_type()){}

    template<unsigned I>
    void apply(){
        result /= closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};

template<unsigned N, class Closure>
typename Closure::value_type abstract_divide(Closure &closure) {
    abstract_subtract_closures<Closure> my_closure(closure);
    meta_loop<N>(my_closure);
    return my_closure.result;
}

///ABS
namespace _my {
template<typename T>
inline T abs(const T& x) {
    return (x < 0) ? -x : x;
}
template<typename T, template<typename Elem, typename = std::allocator<Elem>> class Array = std::vector>
void swap(Array<T>& a, int i, int j) {
    T s = a[i];
    a[i] = a[j];
    a[j] = s;
}

struct print_type
{
    template <class T>
    void operator() (T) const
    {
        auto const& ti = BOOST_CORE_TYPEID(T);
        std::cout << boost::core::demangled_name(ti) << std::endl;
    }
};
template<size_t DimN>
size_t invers(int I, int J, int K) {
    std::array<int, DimN> a{I, J, K};
    size_t count{};
    for(size_t i = 0; i < DimN; ++i) {
        for(size_t j = 0; j < DimN; ++j) {
            if((a[i] > a[j]) && (i < j)) {
                count++;
            }
        }
    }
    return count;
}
template<size_t DimN>
size_t invers(int I, int J, int K, int L) {
    std::array<int, DimN> a{I, J, K, L};
    size_t count{};
    for(size_t i = 1; i < DimN; ++i) {
        for(size_t j = 0; j < DimN; ++j) {
            if((a[i] > a[j]) && (i < j)) {
                count++;
            }
        }
    }
    return count;
}
}


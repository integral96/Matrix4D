#pragma once

#include <include/Matrix3D.hpp>

namespace _spatial {

template<typename Expr>
struct Matrix4D_expr;

struct Matrix4DGrammar : proto::or_<
                              proto::terminal<matrix4_<proto::_>>,
                              proto::plus<Matrix4DGrammar, Matrix4DGrammar>,
                              proto::minus<Matrix4DGrammar, Matrix4DGrammar>,
                              proto::negate< Matrix4DGrammar>,
                              proto::less_equal< Matrix4DGrammar, Matrix4DGrammar>,
                              proto::greater_equal< Matrix4DGrammar, Matrix4DGrammar>,
                              proto::less< Matrix4DGrammar, Matrix4DGrammar>,
                              proto::greater< Matrix4DGrammar, Matrix4DGrammar>,
                              proto::not_equal_to< Matrix4DGrammar, Matrix4DGrammar>,
                              proto::equal_to< Matrix4DGrammar, Matrix4DGrammar>
                                    >{};
struct Matrix4D_domain : proto::domain<proto::generator<Matrix4D_expr>, Matrix4DGrammar> {};

struct Matrix4D_context : proto::callable_context< Matrix4D_context const > {
    Matrix4D_context(size_t i, size_t j, size_t k, size_t l) : i(i), j(j), k(k), l(l) {}

    template<typename Expr, typename Tag = typename Expr::proto_tag>
    struct eval : proto::default_eval<Expr, Matrix4D_context> {};

    template<typename Expr>
    struct eval<Expr, proto::tag::terminal> {
        using result_type = typename proto::result_of::value<Expr>::type::value_type;

        result_type operator()(const Expr& expr, Matrix4D_context& ctx) const {
            return proto::value(expr)(ctx.i, ctx.j, ctx.k, ctx.l);
        }
    };

public:
    size_t i;
    size_t j;
    size_t k;
    size_t l;
};

struct SizeMatrix4D_context {
    SizeMatrix4D_context(size_t Ni, size_t Nj, size_t Nk, size_t Nl)
      : NI(Ni), NJ(Nj), NK(Nk), NL(Nl) {}
    template<typename Expr, typename EnableIf = void>
    struct eval : proto::null_eval<Expr, SizeMatrix4D_context const> {};

    template<typename Expr>
    struct eval<Expr, typename boost::enable_if<
            proto::matches<Expr, proto::terminal<matrix4_<proto::_> > >
        >::type
    >
    {
        typedef void result_type;

        result_type operator ()(Expr &expr, SizeMatrix4D_context const &ctx) const
        {
            if(ctx.NI != proto::value(expr).shape()[0]) {
                throw std::runtime_error("Матрицы не совпадают в размерности по индексу i");
            } else if (ctx.NJ != proto::value(expr).shape()[1]) {
                throw std::runtime_error("Матрицы не совпадают в размерности по индексу j");
            } else if (ctx.NK != proto::value(expr).shape()[2]) {
                throw std::runtime_error("Матрицы не совпадают в размерности по индексу k");
            } else if (ctx.NL != proto::value(expr).shape()[3]) {
                throw std::runtime_error("Матрицы не совпадают в размерности по индексу l");
            }
        }
    };

    size_t NI;
    size_t NJ;
    size_t NK;
    size_t NL;
};

template<typename Expr>
struct Matrix4D_expr : proto::extends<Expr, Matrix4D_expr<Expr>, Matrix4D_domain> {
    Matrix4D_expr(const Expr& expr = Expr()) : Matrix4D_expr::proto_extends(expr) {}

    typename proto::result_of::eval<Expr, Matrix4D_context>::type
    operator () (size_t i, size_t j, size_t k, size_t l) const {
        Matrix4D_context ctx(i, j, k, l);
        return proto::eval(*this, ctx);
    }
};
template<typename T>
struct Matrix4D : Matrix4D_expr<typename proto::terminal< matrix4_<T>>::type> {
    using expr_type = typename proto::terminal< matrix4_<T>>::type;
    using range_tbb = tbb::blocked_rangeNd<size_t, 4>;

    using array_type = typename matrix4_<T>::array_type;
    using range = boost::multi_array_types::index_range;
    using sub_matrix_view1D = typename matrix4_<T>::sub_matrix_view1D;

    const std::array<size_t, 4>& shape_;

    constexpr Matrix4D(const std::array<size_t, 4>& shape) :
        Matrix4D_expr<expr_type>(expr_type::make(matrix4_<T>(shape))), shape_(shape) {

    }
    size_t size(size_t i) const {
        BOOST_ASSERT_MSG((i < 4), "Error i >= 4");
        return proto::value(*this).shape()[i];
    }
    void Random(T min, T max) {
        std::time_t now = std::time(0);
        boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
            if constexpr(std::is_integral_v<T> || IsBigInt<T>::value) {
                boost::random::uniform_int_distribution<> dist{int(min), int(max)};
                for(size_t i = 0; i < size(0); ++i)
                    for(size_t j = 0; j < size(1); ++j)
                        for(size_t k = 0; k < size(2); ++k)
                            for(size_t l = 0; l < size(3); ++l)
                                proto::value(*this)(i, j, k, l) = dist(gen);
            }
            else {
                boost::random::uniform_real_distribution<> dist{double(min), double(max)};
                for(size_t i = 0; i < size(0); ++i)
                    for(size_t j = 0; j < size(1); ++j)
                        for(size_t k = 0; k < size(2); ++k)
                            for(size_t l = 0; l < size(3); ++l)
                                proto::value(*this)(i, j, k, l) = dist(gen);
            }
    }
    template< typename Expr >
    Matrix4D<T> &operator = (Expr const & expr) {
        SizeMatrix4D_context const sizes(size(0), size(1), size(2), size(3));
        proto::eval(proto::as_expr<Matrix4D_domain>(expr), sizes);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
        [&](const range_tbb& out){
        auto out_i = out.dim(0);
        auto out_j = out.dim(1);
        auto out_k = out.dim(2);
        auto out_l = out.dim(3);
        for(size_t i = out_i.begin(); i < out_i.end(); ++i)
            for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                    for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                        proto::value(*this)(i, j, k, l) = expr(i, j, k, l);
        });
        return *this;
    }
    template< typename Expr >
    Matrix4D<T> &operator += (Expr const & expr) {
        SizeMatrix4D_context const sizes(size(0), size(1), size(2), size(3));
        proto::eval(proto::as_expr<Matrix4D_domain>(expr), sizes);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
        [&](const range_tbb& out){
        auto out_i = out.dim(0);
        auto out_j = out.dim(1);
        auto out_k = out.dim(2);
        auto out_l = out.dim(3);
        for(size_t i = out_i.begin(); i < out_i.end(); ++i)
            for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                    for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                        proto::value(*this)(i, j, k, l) += expr(i, j, k, l);
        });
        return *this;
    }
    Matrix4D<T> operator + (const T val) const {
        Matrix4D<T> matrix(shape_);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
        [&](const range_tbb& out){
        auto out_i = out.dim(0);
        auto out_j = out.dim(1);
        auto out_k = out.dim(2);
        auto out_l = out.dim(3);
        for(size_t i = out_i.begin(); i < out_i.end(); ++i)
            for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                    for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                        proto::value(matrix)(i, j, k, l) = proto::value(*this)(i, j, k, l) + val;
        });
        return matrix;
    }
    Matrix4D<T> operator * (const T val) const {
        Matrix4D<T> matrix(shape_);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
        [&](const range_tbb& out){
        auto out_i = out.dim(0);
        auto out_j = out.dim(1);
        auto out_k = out.dim(2);
        auto out_l = out.dim(3);
        for(size_t i = out_i.begin(); i < out_i.end(); ++i)
            for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                    for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                        proto::value(matrix)(i, j, k, l) = proto::value(*this)(i, j, k, l) * val;
        });
        return matrix;
    }
    Matrix4D<T> operator / (const T val) const {
        Matrix4D<T> matrix(shape_);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
        [&](const range_tbb& out){
        auto out_i = out.dim(0);
        auto out_j = out.dim(1);
        auto out_k = out.dim(2);
        auto out_l = out.dim(3);
        for(size_t i = out_i.begin(); i < out_i.end(); ++i)
            for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                    for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                        proto::value(matrix)(i, j, k, l) = proto::value(*this)(i, j, k, l) / val;
        });
        return matrix;
    }
    Matrix3D<T> transversal_matrix(char index, int N) {
        BOOST_ASSERT_MSG((index == 'i') || (index == 'j') || (index == 'k') || (index == 'l'), "Не совпадение индексов");
        typename array_type::index_gen indices;
        if (index == 'i') {
            std::array<size_t, 3> shi{ {size(1), size(2), size(3)} };
            Matrix3D<T> tmp(shi);
            tmp = proto::value(*this)[indices[N][range(0, size(1))][range(0, size(2))][range(0, size(3))]];
            return tmp;
        } else if (index == 'j') {
            std::array<size_t, 3> shj{ {size(0), size(2), size(3)} };
            Matrix3D<T> tmp(shj);
            tmp = proto::value(*this)[indices[range(0, size(0))][N][range(0, size(2))][range(0, size(3))]];
            return tmp;
        } else if (index == 'k') {
            std::array<size_t, 3> shk{ {size(0), size(1), size(3)} };
            Matrix3D<T> tmp(shk);
            tmp = proto::value(*this)[indices[range(0, size(0))][range(0, size(1))][N][range(0, size(3))]];
            return tmp;
        } else if (index == 'l') {
            std::array<size_t, 3> shl{ {size(0), size(1), size(2)} };
            Matrix3D<T> tmp(shl);
            tmp = proto::value(*this)[indices[range(0, size(0))][range(0, size(1))][range(0, size(2))][N]];
            return tmp;
        }
    }
    std::vector<T> transversal_vector(char index, int N)  {
        BOOST_ASSERT_MSG((index == 'i') || (index == 'j') || (index == 'k') || (index == 'l'), "Не совпадение индексов");
        typename array_type::index_gen indices;
        std::vector<T> transversal_vector;
        auto predicate([](int Index, int J, int K, int L) -> bool {
            int j; int M = 4; std::array<int, 4> a{Index, J, K, L};
            do {
                j = Index - 1;
                while ((j != -1) && (a[j] >= a[j + 1])) j--;
                if (j == -1) return false;
                int k = Index - 1;
                while (a[j] >= a[k]) k--;
                _my::swap(a, j, k);
                int l = j + 1, r = Index - 1;
                while (l < r)
                _my::swap(a, l++, r--);
            } while (j > M - 1);
            return true;
        });
        if (index == 'i') {
            for (size_t j = 0; j != size(1); ++j) {
                for (size_t k = 0; k != size(2); ++k) {
                    for (size_t l = 0; l != size(3); ++l) {
                        if(predicate(int(N), int(j), int(k), int(l))){
                            transversal_vector.push_back(transversal_matrix('i', N)(j, k, l));
                        }
                    }
                }
            }
        } else if (index == 'j') {
            for (size_t i = 0; i != size(0); ++i) {
                for (size_t k = 0; k != size(2); ++k) {
                    for (size_t l = 0; l != size(3); ++l) {
                        if(predicate(int(N), int(i), int(k), int(l))){
                            transversal_vector.push_back(transversal_matrix('j', N)(i, k, l));
                        }
                    }
                }
            }
        } else if (index == 'k') {
            for (size_t i = 0; i != size(0); ++i) {
                for (size_t j = 0; j != size(1); ++j) {
                    for (size_t l = 0; l != size(3); ++l) {
                        if(predicate(int(N), int(i), int(j), int(l))){
                            transversal_vector.push_back(transversal_matrix('k', N)(i, j, l));
                        }
                    }
                }
            }
        } else if (index == 'l') {
            for (size_t i = 0; i != size(0); ++i) {
                for (size_t j = 0; j != size(1); ++j) {
                    for (size_t k = 0; k != size(2); ++k) {
                        if(predicate(int(N), int(i), int(j), int(k))){
                            transversal_vector.push_back(transversal_matrix('k', N)(i, j, k));
                        }
                    }
                }
            }
        }
        return transversal_vector;
    }

    T DET_orient(char index, int N) {
        BOOST_ASSERT_MSG((index == 'i') || (index == 'j') || (index == 'k') || (index == 'l'), "Не совпадение индексов");
        if(index == 'i') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(0)), T(1),
                    [=](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector('i', N)[i];
                        } return tmp; }, std::multiplies<T>());
        } else if(index == 'j') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(1)), T(1),
                    [=](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector('j', N)[i];
                        } return tmp; }, std::multiplies<T>());
        } else if(index == 'k') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(2)), T(1),
                    [=](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector('k', N)[i];
                        } return tmp; }, std::multiplies<T>());
        } else if(index == 'l') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(3)), T(1),
                    [=](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector('l', N)[i];
                        } return tmp; }, std::multiplies<T>());
        }
    }
    T DET_FULL() {
        return tbb::parallel_reduce(range_tbb({ 1, size(0) }, { 1, size(1) }, { 1, size(2) }, { 1, size(2) }), T(0),
                [=](const range_tbb& out, T tmp) {
                const auto& out_i = out.dim(0);
                const auto& out_j = out.dim(1);
                const auto& out_k = out.dim(2);
                const auto& out_l = out.dim(3);
                for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                    for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                        for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                            for (size_t l = out_l.begin(); l < out_l.end(); ++l)
                            tmp += std::pow(-1, _my::invers<4>(i, j, k, l))*DET_orient('i', i)*DET_orient('j', j)*DET_orient('k', k)*DET_orient('l', l);
                return tmp; }, std::plus<T>() );
    }

    friend std::ostream& operator << (std::ostream& os, const Matrix4D<T>& A){
            for(const auto& x : proto::value(A)) {
                for(const auto& y : x) {
                    for(const auto& z : y) {
                        for(const auto& f : z) os << f << "\t";
                        os << std::endl;
                    } os << std::endl;
                } os << std::endl;
            } os << std::endl;
            return os;
        }
private:

};

}

#pragma once

#include <include/Base_Matrix.hpp>

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
    using index_type = typename  array_type::index;
    using sub_matrix_view1D = typename matrix4_<T>::sub_matrix_view1D;

    const std::array<index_type, 4>& shape_;

    constexpr Matrix4D(const std::array<index_type, 4>& shape) :
        Matrix4D_expr<expr_type>(expr_type::make(matrix4_<T>(shape))), shape_(shape) {

    }
    size_t size(size_t i) const {
        BOOST_ASSERT_MSG((i < 4), "Error i >= 4");
        return proto::value(*this).shape()[i];
    }
    void Random(T min, T max) {
        std::time_t now = std::time(0);
        boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
            if constexpr(std::is_integral_v<T>) {
                tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
                [&](const range_tbb& out){
                auto out_i = out.dim(0);
                auto out_j = out.dim(1);
                auto out_k = out.dim(2);
                auto out_l = out.dim(3);
                boost::random::uniform_int_distribution<> dist{min, max};
                for(size_t i = out_i.begin(); i < out_i.end(); ++i)
                    for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                        for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                            for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                                proto::value(*this)(i, j, k, l) = dist(gen);
                }, tbb::static_partitioner());
            }
            if constexpr(!std::is_integral_v<T>) {
                tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
                [&](const range_tbb& out){
                auto out_i = out.dim(0);
                auto out_j = out.dim(1);
                auto out_k = out.dim(2);
                auto out_l = out.dim(3);
                boost::random::uniform_real_distribution<> dist{min, max};
                for(size_t i = out_i.begin(); i < out_i.end(); ++i)
                    for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                        for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                            for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                                proto::value(*this)(i, j, k, l) = dist(gen);
                }, tbb::static_partitioner());
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
    decltype(auto) transversal_section(char index) {
        assert ((index == 'i' )|| (index == 'j') || (index == 'k') || (index == 'l'));
        typename array_type::index_gen indices;
        std::vector<sub_matrix_view1D> transversal_vector;
        std::vector<std::array<size_t, 4> > index_vector;
        if(index == 'i'){
            for (index_type j = 0; j != size(1); ++j) {
                for (index_type k = 0; k != size(2); ++k) {
                    for (index_type l = 0; l != size(3); ++l) {
                        if((j + k) == size(0)) {
                            transversal_vector.push_back({proto::value(*this)[ indices[range(0, size(0))][j][k][l] ]});
                            transversal_vector.push_back({proto::value(*this)[ indices[range(0, size(0))][k][j][l] ]});
                            for (index_type i = 0; i != size(0); ++i) {
                                index_vector.push_back({i, j, k, l});
                            }
                            for (index_type i = 0; i != size(0); ++i) {
                                index_vector.push_back({i, k, j, l});
                            }
                        } else if(j == k && (j + k) < 1/*(size_(1) - k)*/) {
                            transversal_vector.push_back({proto::value(*this)[ indices[range(0, size(0))][j][k][l] ]});
                            for (index_type i = 0; i != size(0); ++i) {
                                index_vector.push_back({i, j, k, l});
                            }
                        }
                    }
                }
            }
        }
        if(index == 'j'){
            for (index_type i = 0; i != size(0); ++i) {
                for (index_type k = 0; k != size(2); ++k) {
                    for (index_type l = 0; l != size(2); ++l) {
                        if((i + k) == size(0)) {
                            transversal_vector.push_back({proto::value(*this)[ indices[i][range(0, size(1))][k][l] ]});
                            transversal_vector.push_back({proto::value(*this)[ indices[k][range(0, size(1))][i][l] ]});
                            for (index_type j = 0; j != size(1); ++j) {
                                index_vector.push_back({i, j, k, l});
                            }
                            for (index_type j = 0; j != size(1); ++j) {
                                index_vector.push_back({k, j, i, l});
                            }
                        } else if(i == k && (i + k) < 1/*(size_(1) - k)*/) {
                            transversal_vector.push_back({proto::value(*this)[ indices[i][range(0, size(1))][k][l] ]});
                            for (index_type j = 0; j != size(1); ++j) {
                                index_vector.push_back({i, j, k, l});
                            }
                        }
                    }
                }
            }
        }
        if(index == 'k'){
            for (index_type i = 0; i != size(0); ++i) {
                for (index_type j = 0; j != size(1); ++j) {
                    for (index_type l = 0; l != size(3); ++l) {
                        if((i + j) == size(0)) {
                            transversal_vector.push_back({proto::value(*this)[ indices[i][j][range(0, size(2))][l] ]});
                            transversal_vector.push_back({proto::value(*this)[ indices[j][i][range(0, size(2))][l] ]});
                            for (index_type k = 0; k != size(2); ++k) {
                                index_vector.push_back({i, j, k, l});
                            }
                            for (index_type k = 0; k != size(2); ++k) {
                                index_vector.push_back({j, i, k, l});
                            }
                        } else if(i == j && (i + j) < 1/*(size_(1) - k)*/) {
                            transversal_vector.push_back({proto::value(*this)[ indices[i][j][range(0, size(2))][l] ]});
                            for (index_type k = 0; k != size(2); ++k) {
                                index_vector.push_back({i, j, k, l});
                            }
                        }
                    }
                }
            }
        }
       return std::make_pair(transversal_vector, index_vector);
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

#pragma once

#include <include/Matrix1D.hpp>

namespace _spatial {

template<typename Expr>
struct Matrix2D_expr;

struct Matrix2DGrammar : proto::or_<
                              proto::terminal<matrix2_<proto::_>>,
                              proto::plus<Matrix2DGrammar, Matrix2DGrammar>,
                              proto::minus<Matrix2DGrammar, Matrix2DGrammar>,
                              proto::negate< Matrix2DGrammar>,
                              proto::less_equal< Matrix2DGrammar, Matrix2DGrammar>,
                              proto::greater_equal< Matrix2DGrammar, Matrix2DGrammar>,
                              proto::less< Matrix2DGrammar, Matrix2DGrammar>,
                              proto::greater< Matrix2DGrammar, Matrix2DGrammar>,
                              proto::not_equal_to< Matrix2DGrammar, Matrix2DGrammar>,
                              proto::equal_to< Matrix2DGrammar, Matrix2DGrammar>
                                    >{};
struct Matrix2D_domain : proto::domain<proto::generator<Matrix2D_expr>, Matrix2DGrammar> {};

struct Matrix2D_context : proto::callable_context< Matrix2D_context const > {
    Matrix2D_context(size_t i, size_t j) : i(i), j(j) {}

    template<typename Expr, typename Tag = typename Expr::proto_tag>
    struct eval : proto::default_eval<Expr, Matrix2D_context> {};

    template<typename Expr>
    struct eval<Expr, proto::tag::terminal> {
        using result_type = typename proto::result_of::value<Expr>::type::value_type;

        result_type operator()(const Expr& expr, Matrix2D_context& ctx) const {
            return proto::value(expr)(ctx.i, ctx.j);
        }
    };

public:
    size_t i;
    size_t j;
};

struct SizeMatrix2D_context {
    SizeMatrix2D_context(size_t Ni, size_t Nj)
      : NI(Ni), NJ(Nj) {}
    template<typename Expr, typename EnableIf = void>
    struct eval : proto::null_eval<Expr, SizeMatrix2D_context const> {};

    template<typename Expr>
    struct eval<Expr, typename boost::enable_if<
            proto::matches<Expr, proto::terminal<matrix2_<proto::_> > >
        >::type
    >
    {
        typedef void result_type;

        result_type operator ()(Expr &expr, SizeMatrix2D_context const &ctx) const
        {
            if(ctx.NI != proto::value(expr).shape()[0]) {
                throw std::runtime_error("Матрицы не совпадают в размерности по индексу i");
            } else if (ctx.NJ != proto::value(expr).shape()[1]) {
                throw std::runtime_error("Матрицы не совпадают в размерности по индексу j");
            }
        }
    };

    size_t NI;
    size_t NJ;
};

template<typename Expr>
struct Matrix2D_expr : proto::extends<Expr, Matrix2D_expr<Expr>, Matrix2D_domain> {
    Matrix2D_expr(const Expr& expr = Expr()) : Matrix2D_expr::proto_extends(expr) {}

    typename proto::result_of::eval<Expr, Matrix2D_context>::type
    operator () (size_t i, size_t j) const {
        Matrix2D_context ctx(i, j);
        return proto::eval(*this, ctx);
    }
};
template<typename T>
struct Matrix2D : Matrix2D_expr<typename proto::terminal< matrix2_<T>>::type> {
    using expr_type = typename proto::terminal< matrix2_<T>>::type;
    using range_tbb = tbb::blocked_rangeNd<size_t, 2>;
    using range3_tbb = tbb::blocked_rangeNd<size_t, 3>;

    using array_type = typename matrix2_<T>::array_type;
    using range = boost::multi_array_types::index_range;
    using index_type = typename  array_type::index;
    using sub_matrix_view1D = typename matrix2_<T>::sub_matrix_view1D;

    const std::array<size_t, 2>& shape_;

    constexpr Matrix2D(const std::array<size_t, 2>& shape) :
        Matrix2D_expr<expr_type>(expr_type::make(matrix2_<T>(shape))), shape_(shape) {

    }
    size_t size(size_t i) const {
        BOOST_ASSERT_MSG((i < 2), "Error i >= 2");
        return proto::value(*this).shape()[i];
    }
    void Random(T min, T max) {
        std::time_t now = std::time(0);
        boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
            if constexpr(std::is_integral_v<T>) {
                tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}),
                [&](const range_tbb& out){
                auto out_i = out.dim(0);
                auto out_j = out.dim(1);
                boost::random::uniform_int_distribution<> dist{min, max};
                for(size_t i = out_i.begin(); i < out_i.end(); ++i)
                    for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                                proto::value(*this)(i, j) = dist(gen);
                }, tbb::static_partitioner());
            }
            if constexpr(!std::is_integral_v<T>) {
                tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}),
                [&](const range_tbb& out){
                auto out_i = out.dim(0);
                auto out_j = out.dim(1);
                boost::random::uniform_real_distribution<> dist{min, max};
                for(size_t i = out_i.begin(); i < out_i.end(); ++i)
                    for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                                proto::value(*this)(i, j) = dist(gen);
                }, tbb::static_partitioner());
            }
    }
    template< typename Expr >
    Matrix2D<T> &operator = (Expr const & expr) {
        SizeMatrix2D_context const sizes(size(0), size(1));
        proto::eval(proto::as_expr<Matrix2D_domain>(expr), sizes);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}),
        [&](const range_tbb& out){
        auto out_i = out.dim(0);
        auto out_j = out.dim(1);
        for(size_t i = out_i.begin(); i < out_i.end(); ++i)
            for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                        proto::value(*this)(i, j) = expr(i, j);
        });
        return *this;
    }
    template< typename Expr >
    Matrix2D<T> &operator += (Expr const & expr) {
        SizeMatrix2D_context const sizes(size(0), size(1));
        proto::eval(proto::as_expr<Matrix2D_domain>(expr), sizes);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}),
        [&](const range_tbb& out){
        auto out_i = out.dim(0);
        auto out_j = out.dim(1);
        for(size_t i = out_i.begin(); i < out_i.end(); ++i)
            for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                        proto::value(*this)(i, j) += expr(i, j);
        });
        return *this;
    }
    Matrix2D<T> operator * (Matrix2D<T> const & matr1) {
        SizeMatrix2D_context const sizes(size(0), size(1));
        proto::eval(proto::as_expr<Matrix2D_domain>(matr1), sizes);
        Matrix2D<T> matrix(shape_);
        tbb::parallel_for(range3_tbb({0, size(0)}, {0, size(1)}, {0, size(1)}),
        [&](const range3_tbb& out){
        auto out_i = out.dim(0);
        auto out_j = out.dim(1);
        auto out_k = out.dim(2);
        for(size_t i = out_i.begin(); i < out_i.end(); ++i)
            for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                        proto::value(matrix)(i, j) += proto::value(*this)(i, k)*matr1(k, j);
        });
        return matrix;
    }

    Matrix2D<T> operator + (const T val) const {
        Matrix2D<T> matrix(shape_);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}),
        [&](const range_tbb& out){
        auto out_i = out.dim(0);
        auto out_j = out.dim(1);
        for(size_t i = out_i.begin(); i < out_i.end(); ++i)
            for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                        proto::value(matrix)(i, j) = proto::value(*this)(i, j) + val;
        });
        return matrix;
    }
    Matrix2D<T> operator * (const T val) const {
        Matrix2D<T> matrix(shape_);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}),
        [&](const range_tbb& out){
        auto out_i = out.dim(0);
        auto out_j = out.dim(1);
        for(size_t i = out_i.begin(); i < out_i.end(); ++i)
            for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                        proto::value(matrix)(i, j) = proto::value(*this)(i, j) * val;
        });
        return matrix;
    }

    Matrix2D<T> operator / (const T val) const {
        Matrix2D<T> matrix(shape_);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}),
        [&](const range_tbb& out){
        auto out_i = out.dim(0);
        auto out_j = out.dim(1);
        for(size_t i = out_i.begin(); i < out_i.end(); ++i)
            for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                        proto::value(matrix)(i, j) = proto::value(*this)(i, j) / val;
        });
        return matrix;
    }


    friend std::ostream& operator << (std::ostream& os, const Matrix2D<T>& A){
            for(const auto& x : proto::value(A)) {
                for(const auto& f : x) os << f << "\t";
                os << std::endl;
            } os << std::endl;
            return os;
        }
private:

};

}

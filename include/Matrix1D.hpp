#pragma once

#include "Base_Matrix.hpp"

namespace _spatial {

    template<typename Expr>
    struct Matrix1D_expr;

    struct Matrix1DGrammar : proto::or_<
        proto::terminal<matrix1_<proto::_>>,
        proto::plus<Matrix1DGrammar, Matrix1DGrammar>,
        proto::minus<Matrix1DGrammar, Matrix1DGrammar>,
        proto::negate< Matrix1DGrammar>,
        proto::less_equal< Matrix1DGrammar, Matrix1DGrammar>,
        proto::greater_equal< Matrix1DGrammar, Matrix1DGrammar>,
        proto::less< Matrix1DGrammar, Matrix1DGrammar>,
        proto::greater< Matrix1DGrammar, Matrix1DGrammar>,
        proto::not_equal_to< Matrix1DGrammar, Matrix1DGrammar>,
        proto::equal_to< Matrix1DGrammar, Matrix1DGrammar>
    > {};
    struct Matrix1D_domain : proto::domain<proto::generator<Matrix1D_expr>, Matrix1DGrammar> {};

    struct Matrix1D_context : proto::callable_context< Matrix1D_context const > {
        Matrix1D_context(size_t i) : i(i) {}

        template<typename Expr, typename Tag = typename Expr::proto_tag>
        struct eval : proto::default_eval<Expr, Matrix1D_context> {};

        template<typename Expr>
        struct eval<Expr, proto::tag::terminal> {
            using result_type = typename proto::result_of::value<Expr>::type::value_type;

            result_type operator()(const Expr& expr, Matrix1D_context& ctx) const {
                return proto::value(expr)(ctx.i);
            }
        };

    public:
        size_t i;
    };

    struct SizeMatrix1D_context {
        SizeMatrix1D_context(size_t Ni)
            : NI(Ni) {}
        template<typename Expr, typename EnableIf = void>
        struct eval : proto::null_eval<Expr, SizeMatrix1D_context const> {};

        template<typename Expr>
        struct eval<Expr, typename boost::enable_if<
            proto::matches<Expr, proto::terminal<matrix1_<proto::_> > >
        >::type
        >
        {
            typedef void result_type;

            result_type operator ()(Expr& expr, SizeMatrix1D_context const& ctx) const
            {
                if (ctx.NI != proto::value(expr).shape()[0]) {
                    throw std::runtime_error("Матрицы не совпадают в размерности по индексу i");
                }
            }
        };

        size_t NI;
    };

    template<typename Expr>
    struct Matrix1D_expr : proto::extends<Expr, Matrix1D_expr<Expr>, Matrix1D_domain> {
        Matrix1D_expr(const Expr& expr = Expr()) : Matrix1D_expr::proto_extends(expr) {}

        typename proto::result_of::eval<Expr, Matrix1D_context>::type
            operator () (size_t i) const {
            Matrix1D_context ctx(i);
            return proto::eval(*this, ctx);
        }
    };
    template<typename T>
    struct Matrix1D : Matrix1D_expr<typename proto::terminal< matrix1_<T>>::type> {
        using expr_type = typename proto::terminal< matrix1_<T>>::type;
        using range_tbb = tbb::blocked_rangeNd<size_t, 1>;

        using array_type = typename matrix1_<T>::array_type;
        using range = boost::multi_array_types::index_range;
        using sub_matrix_view1D = typename matrix1_<T>::sub_matrix_view1D;

        const std::array<size_t, 1>& shape_;

        constexpr Matrix1D(const std::array<size_t, 1>& shape) :
            Matrix1D_expr<expr_type>(expr_type::make(matrix1_<T>(shape))), shape_(shape) {

        }
        virtual size_t size() const {
            return proto::value(*this).shape()[0];
        }
        virtual void Random(T min, T max) {
            std::time_t now = std::time(0);
            boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };
            if constexpr (std::is_integral_v<T>) {
                tbb::parallel_for(range_tbb({ 0, size() }),
                    [&](const range_tbb& out) {
                        auto out_i = out.dim(0);
                        boost::random::uniform_int_distribution<> dist{ min, max };
                        for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                            proto::value(*this)(i) = dist(gen);
                    }, tbb::static_partitioner());
            }
            if constexpr (!std::is_integral_v<T>) {
                tbb::parallel_for(range_tbb({ 0, size() }),
                    [&](const range_tbb& out) {
                        auto out_i = out.dim(0);
                        boost::random::uniform_real_distribution<> dist{ min, max };
                        for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                            proto::value(*this)(i) = dist(gen);
                    }, tbb::static_partitioner());
            }
        }
        template< typename Expr >
        Matrix1D<T>& operator = (Expr const& expr) {
            SizeMatrix1D_context const sizes(size());
            proto::eval(proto::as_expr<Matrix1D_domain>(expr), sizes);
            tbb::parallel_for(range_tbb({ 0, size() }),
                [&](const range_tbb& out) {
                    auto out_i = out.dim(0);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        proto::value(*this)(i) = expr(i);
                });
            return *this;
        }
        template< typename Expr >
        Matrix1D<T>& operator += (Expr const& expr) {
            SizeMatrix1D_context const sizes(size());
            proto::eval(proto::as_expr<Matrix1D_domain>(expr), sizes);
            tbb::parallel_for(range_tbb({ 0, size() }),
                [&](const range_tbb& out) {
                    auto out_i = out.dim(0);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        proto::value(*this)(i) += expr(i);
                });
            return *this;
        }
        Matrix1D<T> operator + (const T val) const {
            Matrix1D<T> matrix(shape_);
            tbb::parallel_for(range_tbb({ 0, size() }),
                [&](const range_tbb& out) {
                    auto out_i = out.dim(0);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        proto::value(matrix)(i) = proto::value(*this)(i) + val;
                });
            return matrix;
        }
        Matrix1D<T> operator * (const T val) const {
            Matrix1D<T> matrix(shape_);
            tbb::parallel_for(range_tbb({ 0, size() }),
                [&](const range_tbb& out) {
                    auto out_i = out.dim(0);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        proto::value(matrix)(i) = proto::value(*this)(i) * val;
                });
            return matrix;
        }
        Matrix1D<T> operator / (const T val) const {
            Matrix1D<T> matrix(shape_);
            tbb::parallel_for(range_tbb({ 0, size() }),
                [&](const range_tbb& out) {
                    auto out_i = out.dim(0);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        proto::value(matrix)(i) = proto::value(*this)(i) / val;
                });
            return matrix;
        }


        friend std::ostream& operator << (std::ostream& os, const Matrix1D<T>& A) {
            for (const auto& x : proto::value(A)) {
                os << x << "\t";
            } os << std::endl;
            return os;
        }

    };

}

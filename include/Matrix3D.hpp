#pragma once
#include <include/Matrix2D.hpp>

namespace _spatial {

    template<typename Expr>
    struct Matrix3D_expr;

    struct Matrix3DGrammar : proto::or_<
        proto::terminal<matrix3_<proto::_>>,
        proto::plus<Matrix3DGrammar, Matrix3DGrammar>,
        proto::minus<Matrix3DGrammar, Matrix3DGrammar>,
        proto::negate< Matrix3DGrammar>,
        proto::less_equal< Matrix3DGrammar, Matrix3DGrammar>,
        proto::greater_equal< Matrix3DGrammar, Matrix3DGrammar>,
        proto::less< Matrix3DGrammar, Matrix3DGrammar>,
        proto::greater< Matrix3DGrammar, Matrix3DGrammar>,
        proto::not_equal_to< Matrix3DGrammar, Matrix3DGrammar>,
        proto::equal_to< Matrix3DGrammar, Matrix3DGrammar>
    > {};
    struct Matrix3D_domain : proto::domain<proto::generator<Matrix3D_expr>, Matrix3DGrammar> {};

    struct Matrix3D_context : proto::callable_context< Matrix3D_context const > {
        Matrix3D_context(size_t i, size_t j, size_t k) : i(i), j(j), k(k) {}

        template<typename Expr, typename Tag = typename Expr::proto_tag>
        struct eval : proto::default_eval<Expr, Matrix3D_context> {};

        template<typename Expr>
        struct eval<Expr, proto::tag::terminal> {
            using result_type = typename proto::result_of::value<Expr>::type::value_type;

            result_type operator()(const Expr& expr, Matrix3D_context& ctx) const {
                return proto::value(expr)(ctx.i, ctx.j, ctx.k);
            }
        };

    public:
        size_t i;
        size_t j;
        size_t k;
    };

    struct SizeMatrix3D_context {
        SizeMatrix3D_context(size_t Ni, size_t Nj, size_t Nk)
            : NI(Ni), NJ(Nj), NK(Nk) {}
        template<typename Expr, typename EnableIf = void>
        struct eval : proto::null_eval<Expr, SizeMatrix3D_context const> {};

        template<typename Expr>
        struct eval<Expr, typename boost::enable_if<
            proto::matches<Expr, proto::terminal<matrix3_<proto::_> > >
        >::type
        >
        {
            typedef void result_type;

            result_type operator ()(Expr& expr, SizeMatrix3D_context const& ctx) const
            {
                if (ctx.NI != proto::value(expr).shape()[0]) {
                    throw std::runtime_error("Матрицы не совпадают в размерности по индексу i");
                }
                else if (ctx.NJ != proto::value(expr).shape()[1]) {
                    throw std::runtime_error("Матрицы не совпадают в размерности по индексу j");
                }
                else if (ctx.NK != proto::value(expr).shape()[2]) {
                    throw std::runtime_error("Матрицы не совпадают в размерности по индексу k");
                }
            }
        };

        size_t NI;
        size_t NJ;
        size_t NK;
    };

    template<typename Expr>
    struct Matrix3D_expr : proto::extends<Expr, Matrix3D_expr<Expr>, Matrix3D_domain> {
        Matrix3D_expr(const Expr& expr = Expr()) : Matrix3D_expr::proto_extends(expr) {}

        typename proto::result_of::eval<Expr, Matrix3D_context>::type
            operator () (size_t i, size_t j, size_t k) const {
            Matrix3D_context ctx(i, j, k);
            return proto::eval(*this, ctx);
        }
    };
    template<typename T>
    struct Matrix3D : Matrix3D_expr<typename proto::terminal< matrix3_<T>>::type> {
        using expr_type = typename proto::terminal< matrix3_<T>>::type;
        using range_tbb = tbb::blocked_rangeNd<size_t, 3>;

        using array_type = typename matrix3_<T>::array_type;
        using range = boost::multi_array_types::index_range;
        using sub_matrix_view2D = typename matrix3_<T>::sub_matrix_view2D;
        using sub_matrix_view3D = typename matrix3_<T>::sub_matrix_view3D;

        const std::array<size_t, 3>& shape_;

        constexpr Matrix3D(const std::array<size_t, 3>& shape) :
            Matrix3D_expr<expr_type>(expr_type::make(matrix3_<T>(shape))), shape_(shape) {

        }
        size_t size(size_t i) const {
            BOOST_ASSERT_MSG((i < 3), "Error i >= 4");
            return proto::value(*this).shape()[i];
        }
        void Random(T min, T max) {
            std::time_t now = std::time(0);
            boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
                if constexpr(std::is_integral_v<T>) {
                    boost::random::uniform_int_distribution<> dist{min, max};
                    for(size_t i = 0; i < size(0); ++i)
                        for(size_t j = 0; j < size(1); ++j)
                            for(size_t k = 0; k < size(2); ++k)
                                    proto::value(*this)(i, j, k) = dist(gen);
                }
                if constexpr(!std::is_integral_v<T>) {
                    boost::random::uniform_real_distribution<> dist{min, max};
                    for(size_t i = 0; i < size(0); ++i)
                        for(size_t j = 0; j < size(1); ++j)
                            for(size_t k = 0; k < size(2); ++k)
                                    proto::value(*this)(i, j, k) = dist(gen);
                }
        }
        template< typename Expr >
        Matrix3D<T>& operator = (Expr const& expr) {
            SizeMatrix3D_context const sizes(size(0), size(1), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(expr), sizes);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(*this)(i, j, k) = expr(i, j, k);
                });
            return *this;
        }
        Matrix3D<T>& operator = (sub_matrix_view3D const& expr) {
            SizeMatrix3D_context const sizes(size(0), size(1), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(expr), sizes);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(*this)(i, j, k) = expr[i][j][k];
                });
            return *this;
        }
        template< typename Expr >
        Matrix3D<T>& operator += (Expr const& expr) {
            SizeMatrix3D_context const sizes(size(0), size(1), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(expr), sizes);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(*this)(i, j, k) += expr(i, j, k);
                });
            return *this;
        }
        Matrix3D<T> operator + (const T val) const {
            Matrix3D<T> matrix(shape_);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(matrix)(i, j, k) = proto::value(*this)(i, j, k) + val;
                });
            return matrix;
        }
        Matrix3D<T> operator * (const T val) const {
            Matrix3D<T> matrix(shape_);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(matrix)(i, j, k) = proto::value(*this)(i, j, k) * val;
                });
            return matrix;
        }
        Matrix3D<T> operator / (const T val) const {
            Matrix3D<T> matrix(shape_);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(matrix)(i, j, k) = proto::value(*this)(i, j, k) / val;
                });
            return matrix;
        }

        std::vector<Matrix2D<T>> transversal(char index)  {
            BOOST_ASSERT_MSG((index == 'i') || (index == 'j') || (index == 'k') , "Не совпадение индексов");
            typename array_type::index_gen indices;
            std::array<size_t, 2> shi{ {size(1), size(2)} };
            std::array<size_t, 2> shj{ {size(0), size(2)} };
            std::array<size_t, 2> shk{ {size(0), size(1)} };
            std::vector<Matrix2D<T>> transversal_vector;
            if (index == 'i') {
                for (size_t i = 0; i != size(0); ++i) {
                    auto tmp = std::make_unique<Matrix2D<T>>(shi);
                    *tmp = proto::value(*this)[indices[i][range(0, size(1))][range(0, size(2))]];
                    transversal_vector.push_back(*tmp);
                }
            }
            else if (index == 'j') {
                for (size_t j = 0; j != size(1); ++j) {
                    auto tmp = std::make_unique<Matrix2D<T>>(shj);
                    *tmp = proto::value(*this)[indices[range(0, size(0))][j][range(0, size(2))]];
                    transversal_vector.push_back(*tmp);
                }
            }
            else if (index == 'k') {
                for (size_t k = 0; k != size(2); ++k) {
                    auto tmp = std::make_unique<Matrix2D<T>>(shk);
                    *tmp = proto::value(*this)[indices[range(0, size(0))][range(0, size(1))][k]];
                    transversal_vector.push_back(*tmp);
                }
            }
            return transversal_vector;
        }
        int64_t DET_orient(char index) {
            BOOST_ASSERT_MSG((index == 'i') || (index == 'j') || (index == 'k') , "Не совпадение индексов");
            if(index == 'i') {
                return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(0)), 1.0,
                        [=](const tbb::blocked_range<size_t>& r, T tmp) {
                            for (size_t i = r.begin(); i != r.end(); ++i) {
                                tmp *= transversal('i')[i].DET();
                            } return tmp; }, std::multiplies<T>());
            } else if(index == 'j') {
                return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(1)), 1.0,
                        [=](const tbb::blocked_range<size_t>& r, T tmp) {
                            for (size_t i = r.begin(); i != r.end(); ++i) {
                                tmp *= transversal('j')[i].DET();
                            } return tmp; }, std::multiplies<T>());
            } else if(index == 'k') {
                return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(2)), 1.0,
                        [=](const tbb::blocked_range<size_t>& r, T tmp) {
                            for (size_t i = r.begin(); i != r.end(); ++i) {
                                tmp *= transversal('k')[i].DET();
                            } return tmp; }, std::multiplies<T>());
            }
        }
        int64_t DET_FULL() {
            return DET_orient('i')*DET_orient('j')*DET_orient('k');
        }

        friend std::ostream& operator << (std::ostream& os, const Matrix3D<T>& A) {
            for (const auto& x : proto::value(A)) {
                for (const auto& y : x) {
                    for (const auto& f : y) os << f << "\t";
                    os << std::endl;
                } os << std::endl;
            } os << std::endl;
            return os;
        }
    private:

    };

}

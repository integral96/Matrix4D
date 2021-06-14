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

    ///Product struct matrix
    template<class MatrND, class MatrN_1D>
    struct product_closure {
    private:
        const MatrND& A;
        const MatrN_1D& a;
    public:
        product_closure(const MatrND& A_, const MatrN_1D& a_) :
         A(A_), a(a_) {}
        template<char index_name>
        T value(size_t K, size_t i, size_t j, size_t k, size_t l) const {
            if constexpr(IsMatrix3D<T, MatrND>::value) {
                if constexpr(index_name == 'i') {
                    return proto::value(A)(K, j, k)*proto::value(a)(K, i, l);
                }
                if constexpr(index_name == 'j') {
                    return proto::value(A)(i, K, k)*proto::value(a)(K, j, l);
                }
                if constexpr(index_name == 'k') {
                    return proto::value(A)(i, j, K)*proto::value(a)(K, k, l);
                }
            } else if constexpr(IsMatrix4D<T, MatrND>::value) {
                if constexpr(index_name == 'i') {
                    return proto::value(A)(K, j, k, l)*proto::value(a)(K, i, l);
                }
                if constexpr(index_name == 'j') {
                    return proto::value(A)(i, K, k, l)*proto::value(a)(K, j, l);
                }
                if constexpr(index_name == 'k') {
                    return proto::value(A)(i, j, K, l)*proto::value(a)(K, k, l);
                }
                if constexpr(index_name == 'l') {
                    return proto::value(A)(i, j, k, K)*proto::value(a)(K, l, l);
                }
            }
        }
    };
    template<class MatrND, class MatrN_1D>
    struct product_ {
    private:
        Matrix4D<T>& this_;
        const MatrND& A;
        const MatrN_1D& a;
        T result;
        template<char index_name>
        const T summ_(product_closure<MatrND, MatrN_1D>& closure, size_t i, size_t j, size_t k, size_t l) {
            if constexpr(IsMatrix3D<T, MatrND>::value) {
                if constexpr(index_name == 'i') {
                    for(size_t K = 0; K < proto::value(this_).shape()[0]; ++K)
                        result += closure.template value<'i'>(K, i, j, k, l);
                    return result;
                }
                if constexpr(index_name == 'j') {
                    for(size_t K = 0; K < proto::value(this_).shape()[1]; ++K)
                        result += closure.template value<'j'>(K, i, j, k, l);
                    return result;
                }
                if constexpr(index_name == 'k') {
                    for(size_t K = 0; K < proto::value(this_).shape()[2]; ++K)
                        result += closure.template value<'k'>(K, i, j, k, l);
                    return result;
                }
            } else if constexpr(IsMatrix4D<T, MatrND>::value) {
                if constexpr(index_name == 'i') {
                    for(size_t K = 0; K < proto::value(this_).shape()[0]; ++K)
                        result += closure.template value<'i'>(K, i, j, k, l);
                    return result;
                }
                if constexpr(index_name == 'j') {
                    for(size_t K = 0; K < proto::value(this_).shape()[1]; ++K)
                        result += closure.template value<'j'>(K, i, j, k, l);
                    return result;
                }
                if constexpr(index_name == 'k') {
                    for(size_t K = 0; K < proto::value(this_).shape()[2]; ++K)
                        result += closure.template value<'k'>(K, i, j, k, l);
                    return result;
                }
                if constexpr(index_name == 'l') {
                    for(size_t K = 0; K < proto::value(this_).shape()[3]; ++K)
                        result += closure.template value<'l'>(K, i, j, k, l);
                    return result;
                }
            }
        }
    public:
        product_(Matrix4D<T>& this_, const MatrND& A_, const MatrN_1D& a_) :
        this_(this_), A(A_), a(a_) {}
        template<char index_name>
        void init(size_t i, size_t j, size_t k, size_t l) {
            product_closure<MatrND, MatrN_1D> closure(A, a);
            proto::value(this_)(i, j, k, l) = summ_<index_name>(closure, i, j, k, l);
        }
    };

    ///Finish struct product

    const std::array<size_t, 4>& shape_;

    constexpr Matrix4D(const std::array<size_t, 4>& shape) :
        Matrix4D_expr<expr_type>(expr_type::make(matrix4_<T>(shape))), shape_(shape) {

    }
    size_t size(size_t i) const {
        BOOST_ASSERT_MSG((i < 4), "Error i >= 4");
        return proto::value(*this).shape()[i];
    }
    void init(const std::vector<std::vector<std::vector<std::vector<T>>>>& list) {
        BOOST_ASSERT_MSG(list.size() == size(0), "size orient i not equal");
        BOOST_ASSERT_MSG(list.begin()->size() == size(1), "size orient j not equal");
        BOOST_ASSERT_MSG(list.begin()->begin()->size() == size(2), "size orient k not equal");
        BOOST_ASSERT_MSG(list.begin()->begin()->begin()->size() == size(2), "size orient k not equal");
        for(size_t i = 0; i < size(0); ++i)
            for(size_t j = 0; j < size(1); ++j)
                for(size_t k = 0; k < size(2); ++k)
                    for(size_t l = 0; l < size(3); ++l)
                        proto::value(*this)(i, j, k, l) = list.at(i).at(j).at(k).at(l);
    }
    void Random(T min, T max, long shift = 0) {
        std::time_t now = std::time(&shift);
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
    ///Product 3D
    template<char index_name>
    Matrix4D<T>& prod (Matrix3D<T> const& A, Matrix3D<T> const& B) {
        if constexpr(index_name == 'i'){
            SizeMatrix3D_context const sizes(size(0), size(0), size(0));
            proto::eval(proto::as_expr<Matrix3D_domain>(A), sizes);
            proto::eval(proto::as_expr<Matrix3D_domain>(B), sizes);
        }
        if constexpr(index_name == 'j'){
            SizeMatrix3D_context const sizes(size(1), size(1), size(1));
            proto::eval(proto::as_expr<Matrix3D_domain>(A), sizes);
            proto::eval(proto::as_expr<Matrix3D_domain>(B), sizes);
        }
        if constexpr(index_name == 'k'){
            SizeMatrix3D_context const sizes(size(2), size(2), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(A), sizes);
            proto::eval(proto::as_expr<Matrix3D_domain>(B), sizes);
        }
        tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }, { 0, size(3) }),
            [&](const range_tbb& out) {
                const auto& out_i = out.dim(0);
                const auto& out_j = out.dim(1);
                const auto& out_k = out.dim(2);
                const auto& out_l = out.dim(3);
                for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                    for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                        for (size_t k = out_k.begin(); k < out_k.end(); ++k) {
                            for (size_t l = out_l.begin(); l < out_l.end(); ++l) {
                                product_<Matrix3D<T>, Matrix3D<T>>(*this, A, B).template init<index_name>(i, j, k, l);
                            }
                        }
            });
        return *this;
    }
    ///Product 4D
    template<char index_name>
    Matrix4D<T>& prod (Matrix4D<T> const& A, Matrix3D<T> const& B) {
        if constexpr(index_name == 'i'){
            SizeMatrix3D_context const sizes(size(0), size(0), size(0));
            proto::eval(proto::as_expr<Matrix3D_domain>(B), sizes);
        }
        if constexpr(index_name == 'j'){
            SizeMatrix3D_context const sizes(size(1), size(1), size(1));
            proto::eval(proto::as_expr<Matrix3D_domain>(B), sizes);
        }
        if constexpr(index_name == 'k'){
            SizeMatrix3D_context const sizes(size(2), size(2), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(B), sizes);
        }
        SizeMatrix4D_context const sizes(size(0), size(1), size(2), size(3));
        proto::eval(proto::as_expr<Matrix4D_domain>(A), sizes);
        tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }, { 0, size(3) }),
            [&](const range_tbb& out) {
                const auto& out_i = out.dim(0);
                const auto& out_j = out.dim(1);
                const auto& out_k = out.dim(2);
                const auto& out_l = out.dim(3);
                for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                    for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                        for (size_t k = out_k.begin(); k < out_k.end(); ++k) {
                            for (size_t l = out_l.begin(); l < out_l.end(); ++l) {
                                product_<Matrix4D<T>, Matrix3D<T>>(*this, A, B).template init<index_name>(i, j, k, l);
                            }
                        }
            });
        return *this;
    }
    template<char index>
    Matrix3D<T> transversal_matrix(int N) {
        BOOST_STATIC_ASSERT_MSG((index == 'i') || (index == 'j') || (index == 'k') || (index == 'l'), "Не совпадение индексов");
        typename array_type::index_gen indices;
        if constexpr(index == 'i') {
            std::array<size_t, 3> shi{ {size(1), size(2), size(3)} };
            Matrix3D<T> tmp(shi);
            tmp = proto::value(*this)[indices[N][range(0, size(1))][range(0, size(2))][range(0, size(3))]];
            return tmp;
        } else if constexpr(index == 'j') {
            std::array<size_t, 3> shj{ {size(0), size(2), size(3)} };
            Matrix3D<T> tmp(shj);
            tmp = proto::value(*this)[indices[range(0, size(0))][N][range(0, size(2))][range(0, size(3))]];
            return tmp;
        } else if constexpr(index == 'k') {
            std::array<size_t, 3> shk{ {size(0), size(1), size(3)} };
            Matrix3D<T> tmp(shk);
            tmp = proto::value(*this)[indices[range(0, size(0))][range(0, size(1))][N][range(0, size(3))]];
            return tmp;
        } else if constexpr(index == 'l') {
            std::array<size_t, 3> shl{ {size(0), size(1), size(2)} };
            Matrix3D<T> tmp(shl);
            tmp = proto::value(*this)[indices[range(0, size(0))][range(0, size(1))][range(0, size(2))][N]];
            return tmp;
        }
    }
    template<char index>
    std::vector<T> transversal_vector(int N)  {
        BOOST_STATIC_ASSERT_MSG((index == 'i') || (index == 'j') || (index == 'k') || (index == 'l'), "Не совпадение индексов");
        typename array_type::index_gen indices;
        std::vector<T> transversal_vector;
        transversal_vector.reserve(size(0)*size(1)*size(2)*size(3));
        struct index4D {  int I; int J; int K; int L; };
        auto predicate([s_ = std::max({size(0), size(1), size(2), size(3)})]
                       (std::vector<index4D>& index_) {
            for(int i = 0; i < s_; ++i) {
                if((index_.at(i).I != index_.at(i + 1).I) &&
                        (index_.at(i).J == index_.at(i + 1).K ||
                         index_.at(i).K == index_.at(i + 1).J ||
                         index_.at(i).J == index_.at(i + 1).L)) return true;
                else return false;
            }
        });

        if constexpr(index == 'i') {
            std::vector<index4D> index_vector;
            index_vector.reserve(size(0)*size(1)*size(2)*size(3));
            for (size_t j = 0; j != size(1); ++j) {
                for (size_t k = 0; k != size(2); ++k) {
                    for (size_t l = 0; l != size(3); ++l) {
                        index_vector.push_back({(int)N, int(j), int(k), int(l)});
                        index_vector.push_back({(int)(N + 1), int(j), int(k), int(l)});
                        if(predicate(index_vector)){
                            transversal_vector.push_back(transversal_matrix<index>(N)(j, k, l));
                        }
                    }
                }
            }
        } else if constexpr(index == 'j') {
            std::vector<index4D> index_vector;
            index_vector.reserve(size(0)*size(1)*size(2)*size(3));
            for (size_t i = 0; i != size(0); ++i) {
                for (size_t k = 0; k != size(2); ++k) {
                    for (size_t l = 0; l != size(3); ++l) {
                        index_vector.push_back({(int)N, int(i), int(k), int(l)});
                        index_vector.push_back({(int)(N + 1), int(i), int(k), int(l)});
                        if(predicate(index_vector)){
                            transversal_vector.push_back(transversal_matrix<index>(N)(i, k, l));
                        }
                    }
                }
            }
        } else if constexpr(index == 'k') {
            std::vector<index4D> index_vector;
            index_vector.reserve(size(0)*size(1)*size(2)*size(3));
            for (size_t i = 0; i != size(0); ++i) {
                for (size_t j = 0; j != size(1); ++j) {
                    for (size_t l = 0; l != size(3); ++l) {
                        index_vector.push_back({(int)N, int(i), int(j), int(l)});
                        index_vector.push_back({(int)(N + 1), int(i), int(j), int(l)});
                        if(predicate(index_vector)){
                            transversal_vector.push_back(transversal_matrix<index>(N)(i, j, l));
                        }
                    }
                }
            }
        } else if constexpr(index == 'l') {
            std::vector<index4D> index_vector;
            index_vector.reserve(size(0)*size(1)*size(2)*size(3));
            for (size_t i = 0; i != size(0); ++i) {
                for (size_t j = 0; j != size(1); ++j) {
                    for (size_t k = 0; k != size(2); ++k) {
                        index_vector.push_back({(int)N, int(i), int(j), int(k)});
                        index_vector.push_back({(int)(N + 1), int(i), int(j), int(k)});
                        if(predicate(index_vector)){
                            transversal_vector.push_back(transversal_matrix<index>(N)(i, j, k));
                        }
                    }
                }
            }
        }
        return transversal_vector;
    }
    template<char index>
    T DET_orient(int N) {
        BOOST_STATIC_ASSERT_MSG((index == 'i') || (index == 'j') || (index == 'k') || (index == 'l'), "Не совпадение индексов");
        if constexpr(index == 'i') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(0)), T(1),
                    [this, &N](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector<index>(N)[i];
                        } return tmp; }, std::multiplies<T>());
        } else if constexpr(index == 'j') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(1)), T(1),
                    [this, &N](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector<index>(N)[i];
                        } return tmp; }, std::multiplies<T>());
        } else if constexpr(index == 'k') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(2)), T(1),
                    [this, &N](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector<index>(N)[i];
                        } return tmp; }, std::multiplies<T>());
        } else if constexpr(index == 'l') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(3)), T(1),
                    [this, &N](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector<index>(N)[i];
                        } return tmp; }, std::multiplies<T>());
        }
    }
    T DET_FULL() {
        return tbb::parallel_reduce(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }, { 0, size(3) }), T(0),
                [this](const range_tbb& out, T tmp) {
                const auto& out_i = out.dim(0);
                const auto& out_j = out.dim(1);
                const auto& out_k = out.dim(2);
                const auto& out_l = out.dim(3);
                for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                    for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                        for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                            for (size_t l = out_l.begin(); l < out_l.end(); ++l)
                            tmp += std::pow(-1, _my::invers<4>(i, j, k, l))*DET_orient<'i'>(i)*DET_orient<'j'>(j)*DET_orient<'k'>(k)*DET_orient<'l'>(l);
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

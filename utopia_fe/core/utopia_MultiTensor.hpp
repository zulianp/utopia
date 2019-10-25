#ifndef UTOPIA_FE_TENSOR_HPP
#define UTOPIA_FE_TENSOR_HPP

#include "moonolith_expanding_array.hpp"
#include "utopia_Traits.hpp"
#include "utopia_blas.hpp"

#include <iostream>
#include <cmath>

namespace utopia {

    static const int FEBLAS = BLAS + 1;

    template<class T_, int Order_>
    class MultiTensor;

    template<typename T>
    class BlasMultiTensorTraits {
    public:
        using BlasTraitsT = utopia::Traits<BlasVector<T>>;
        using Scalar      = utopia::MultiTensor<T, 0>;
        using ScalarValue = T;
        using Matrix      = utopia::MultiTensor<utopia::BlasMatrix<T>, 2>;
        using Vector      = utopia::MultiTensor<utopia::BlasVector<T>, 1>;

        using SizeType    = typename BlasTraitsT::SizeType;
        using IndexSet    = typename BlasTraitsT::IndexSet;
        using ScalarArray = typename BlasTraitsT::ScalarArray;
        using IndexArray  = typename BlasTraitsT::IndexArray;

        enum {
            Backend = FEBLAS
        };

        static BackendInfo &backend_info()
        {
            static BackendInfo instance_("fe_blas");
            return instance_;
        }
    };

    using BlasMultiTensorTraitsd = utopia::BlasMultiTensorTraits<double>;
    using MultiVectord = utopia::MultiTensor<BlasVectord, 1>;
    using MultiMatrixd = utopia::MultiTensor<BlasMatrixd, 2>;
    using MultiScalard = utopia::MultiTensor<double, 0>;

    UTOPIA_MAKE_TRAITS(MultiScalard, BlasMultiTensorTraitsd, 0);
    UTOPIA_MAKE_TRAITS(MultiVectord, BlasMultiTensorTraitsd, 1);
    UTOPIA_MAKE_TRAITS_DENSE(MultiMatrixd, BlasMultiTensorTraitsd, 2);


    template<typename T>
    class Math<MultiTensor<T, 0>> {
    public:
        inline static MultiTensor<T, 0> abs(const MultiTensor<T, 0> &x) { 
            auto n = x.size();
            MultiTensor<T, 0> ret = x;
            for(SizeType i = 0; i < n; ++i) {
                ret[i] = std::abs(ret[i]);
            }
            
            return ret;
        }

        inline static MultiTensor<T, 0> abs(MultiTensor<T, 0> &&x) { 
            auto n = x.size();
            MultiTensor<T, 0> ret = x;
            for(SizeType i = 0; i < n; ++i) {
                ret[i] = std::abs(ret[i]);
            }
            
            return ret;
        }

    };

    template<typename T>
    MultiTensor<T, 0> operator*(const MultiTensor<T, 0> &left, const MultiTensor<T, 0> &right)
    {
        auto n = left.size();
        MultiTensor<T, 0> ret = left;
        for(SizeType i = 0; i < n; ++i) {
            ret[i] *= right[i];
        }

        return ret;
    }

    template<class T_, int Order_ = Traits<T_>::Order>
    class MultiTensor : public Tensor<MultiTensor<T_>, Order_> {
    public:
        using T = T_;
        static const int Order = Order_;
        
        using SizeType = typename Traits<MultiTensor>::SizeType;
        using Scalar   = typename Traits<MultiTensor>::Scalar;
        using ScalarValue = typename Traits<MultiTensor>::ScalarValue;
        
        
        using Super = utopia::Tensor<MultiTensor, Order_>;
        using Super::Super;

        MultiTensor(MultiTensor &&other) 
        : values_(std::move(other.values_))
        {}

        MultiTensor(const MultiTensor &other) 
        : values_(other.values_)
        {}

        MultiTensor &operator=(MultiTensor &&other) 
        {
            values_ = std::move(other.values_);
            return *this;
        }

        MultiTensor &operator=(const MultiTensor &other) 
        {
            if(is_alias(other)) return *this;
            values_ = other.values_;
            return *this;
        }

        template<class Expr>
        MultiTensor(const Expression<Expr> &expr)
        {
            //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::construct_eval(expr.derived());
        }

        template<class Expr>
        inline MultiTensor &operator=(const Expression<Expr> &expr)
        {
            //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::assign_eval(expr.derived());
            return *this;
        }

        void assign(const MultiTensor &other) override
        {
            values_ = other.values_;
        }

        void assign(MultiTensor &&other) override
        {
            values_ = std::move(other.values_);
        }

        template<class TAlpha, int OrderAlpha>
        void axpy(
            const MultiTensor<TAlpha, OrderAlpha> &alpha,
            const MultiTensor &x)
        {
            const auto n = size();
            assert(n == alpha.size());
            assert(n == x.size());

            for(SizeType i = 0; i < n; ++i) {
                at(i) += alpha[i] * x[i];
            }
        }

        void axpy(const ScalarValue &alpha, const MultiTensor &x)
        {
            const auto n = size();
            assert(n == x.size());

            for(SizeType i = 0; i < n; ++i) {
                at(i) += alpha * x[i];
            }
        }

        void dot(const MultiTensor &other, Scalar &result) const
        {
            const auto n = size();
            assert(n == other.size());

            result.resize(n);

            for(SizeType i = 0; i < n; ++i) {
                result[i] = utopia::dot(at(i), other[i]);
            }
        }

        Scalar dot(const MultiTensor &other) const
        {
            Scalar result;
            dot(other, result);
            return result;
        }

        void scale(const Scalar &alpha)
        {
            const auto n = size();
            assert(n == alpha.size());

            for(SizeType i = 0; i < n; ++i) {
                at(i) *= alpha[i];
            }
        }

        void scale(const ScalarValue &alpha)
        {
            const auto n = size();
            assert(n == alpha.size());

            for(SizeType i = 0; i < n; ++i) {
                at(i) *= alpha;
            }
        }

        inline SizeType size() const
        {
            return values_.size();
        }

        inline void resize(const SizeType n)
        {
            values_.resize(n);
        }

        inline const T &operator[](const SizeType i) const
        {
            assert(i < size());
            return values_[i];
        }

        inline T &operator[](const SizeType i)
        {
            assert(i < size());
            return values_[i];
        }

        inline T &at(const SizeType i)
        {
            assert(i < size());
            return values_[i];
        }

        inline const T &at(const SizeType i) const
        {
            assert(i < size());
            return values_[i];
        }

        inline bool is_alias(const MultiTensor &other) const
        {
            return this == &other;
            //if we move to views this is actually better
            //return &values_[0] == &other.values_[0];
        }

        inline void describe() const
        {
            const auto n = size();

            std::cout << "----------------------\n";
            for(SizeType i = 0; i < n; ++i) {
                std::cout << i << ")\n";
                disp(at(i));
            }
            std::cout << "----------------------\n";
        }

     
    private:
        moonolith::Storage<T> values_;
    };
}

#endif //UTOPIA_FE_TENSOR_HPP

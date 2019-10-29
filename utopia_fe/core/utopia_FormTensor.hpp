#ifndef UTOPIA_FORM_TENSOR_HPP
#define UTOPIA_FORM_TENSOR_HPP

#include "moonolith_expanding_array.hpp"
#include "utopia_LocalFEForwardDeclarations.hpp"
#include "utopia_MultiTensor.hpp"

namespace utopia {

    template<typename T>
    class BlasFormTensorTraits {
    public:
        using BlasTraitsT = utopia::Traits<BlasVector<T>>;
        using Scalar = T;
        using MultiScalar = utopia::FormTensor<T, 0>;
        using Matrix      = utopia::FormTensor<utopia::BlasMatrix<T>, 2>;
        using Vector      = utopia::FormTensor<utopia::BlasVector<T>, 1>;

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

    using BlasFormTensorTraitsd = utopia::BlasFormTensorTraits<double>;
    using FormVectord = utopia::FormTensor<BlasVectord, 1>;
    using FormMatrixd = utopia::FormTensor<BlasMatrixd, 2>;
    using FormScalard = utopia::FormTensor<double, 0>;

    UTOPIA_MAKE_TRAITS(FormScalard, BlasFormTensorTraitsd, 0);
    UTOPIA_MAKE_TRAITS(FormVectord, BlasFormTensorTraitsd, 1);
    UTOPIA_MAKE_TRAITS_DENSE(FormMatrixd, BlasFormTensorTraitsd, 2);


    template<typename T2, int Order>
    Binary<Number<typename FormTensor<T2, Order>::Scalar>, FormTensor<T2, Order>, Multiplies> 
    operator*(const typename FormTensor<T2, Order>::Scalar &left, const FormTensor<T2, Order> &right)
    {
        return Binary<Number<typename FormTensor<T2, Order>::Scalar>, FormTensor<T2, Order>, Multiplies>(left, right);
    }

    template<class T_, int Order_>
    class FormTensor : public Tensor<FormTensor<T_, Order_>, Order_> {
    public:
        using T = T_;
        static const int Order = Order_;

        using SizeType    = typename Traits<FormTensor>::SizeType;
        using MultiScalar = typename Traits<FormTensor>::MultiScalar;
        using Scalar      = typename Traits<FormTensor>::Scalar;
        
        using Super = utopia::Tensor<FormTensor, Order_>;
        using Super::Super;

        inline std::string get_class() const override
        {
            return "Form" + Order2String<Order>::Value();
        }

        FormTensor(FormTensor &&other) 
        : values_(std::move(other.values_))
        {}

        FormTensor(const FormTensor &other) 
        : values_(other.values_)
        {}

        FormTensor &operator=(FormTensor &&other) 
        {
            values_ = std::move(other.values_);
            return *this;
        }

        FormTensor &operator=(const FormTensor &other) 
        {
            if(is_alias(other)) return *this;
            values_ = other.values_;
            return *this;
        }

        template<class Expr>
        FormTensor(const Expression<Expr> &expr)
        {
            //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::construct_eval(expr.derived());
        }

        template<class Expr>
        inline FormTensor &operator=(const Expression<Expr> &expr)
        {
            //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::assign_eval(expr.derived());
            return *this;
        }

        void assign(const FormTensor &other) override
        {
            values_ = other.values_;
        }

        void assign(FormTensor &&other) override
        {
            values_ = std::move(other.values_);
        }

        void resize(const SizeType n)
        {
            values_.resize(n);
        }

        inline SizeType size() const
        {
            return values_.size();
        }

        inline const MultiTensor<T, Order> &operator[](const SizeType i) const
        {
            assert(i < size());
            return values_[i];
        }

        inline MultiTensor<T, Order> &operator[](const SizeType i)
        {
            assert(i < size());
            return values_[i];
        }

        inline MultiTensor<T, Order> &at(const SizeType i)
        {
            assert(i < size());
            return values_[i];
        }

        inline const MultiTensor<T, Order> &at(const SizeType i) const
        {
            assert(i < size());
            return values_[i];
        }

        void describe() const
        {
            SizeType idx = 0;
            for(const auto &v : values_)
            {
                std::cout << idx++ << "]\n";
                disp(v);
            }
        }

        inline void scale(const Scalar &factor)
        {
            for(auto &v : values_) {
                v.scale(factor);
            }
        }

        inline bool is_alias(const FormTensor &other) const
        {
            return this == &other;
        }

    private:
        moonolith::Storage<MultiTensor<T, Order>> values_;
    };
}

#endif //UTOPIA_FORM_TENSOR_HPP
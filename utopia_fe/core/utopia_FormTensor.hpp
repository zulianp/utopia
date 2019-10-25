#ifndef UTOPIA_FORM_TENSOR_HPP
#define UTOPIA_FORM_TENSOR_HPP

#include "moonolith_expanding_array.hpp"
#include "utopia_MultiTensor.hpp"


namespace utopia {

    template<class T_, int Order>
    class FormTensor;

    template<typename T>
    class BlasFormTensorTraits {
    public:
        using BlasTraitsT = utopia::Traits<BlasVector<T>>;
        using Scalar      = utopia::FormTensor<T, 0>;
        using ScalarValue = T;
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

    using BlasFormTensorTraitsd = utopia::BlasMultiTensorTraits<double>;
    using FormVectord = utopia::FormTensor<BlasVectord, 1>;
    using FormMatrixd = utopia::FormTensor<BlasMatrixd, 2>;
    using FormScalard = utopia::FormTensor<double, 0>;

    UTOPIA_MAKE_TRAITS(FormScalard, BlasFormTensorTraitsd, 0);
    UTOPIA_MAKE_TRAITS(FormVectord, BlasFormTensorTraitsd, 1);
    UTOPIA_MAKE_TRAITS_DENSE(FormMatrixd, BlasFormTensorTraitsd, 2);

    template<class T_, int Order_>
    class FormTensor : public Tensor<FormTensor<T_, Order_>, Order_> {
    public:
        using T = T_;
        static const int Order = Order_;

        using SizeType = typename Traits<FormTensor>::SizeType;
        using Scalar   = typename Traits<FormTensor>::Scalar;
        using ScalarValue = typename Traits<FormTensor>::ScalarValue;
        
        
        using Super = utopia::Tensor<FormTensor, Order_>;
        using Super::Super;

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

    private:
        moonolith::Storage<MultiTensor<T, Order>> values_;
    };
}

#endif //UTOPIA_FORM_TENSOR_HPP
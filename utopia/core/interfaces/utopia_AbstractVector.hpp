#ifndef UTOPIA_ABSTRACT_VECTOR_HPP
#define UTOPIA_ABSTRACT_VECTOR_HPP

#include "utopia_Base.hpp"
#include "utopia_Range.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_Vector.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_Normed.hpp"
#include "utopia_Reducible.hpp"
#include "utopia_Transformable.hpp"
#include "utopia_ElementWiseOperand.hpp"
#include "utopia_Constructible.hpp"
#include "utopia_Comparable.hpp"
#include "utopia_BLAS_Operands.hpp"
#include "utopia_Allocations.hpp"
#include "utopia_Select.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {

    //parallel types, collective operations
    template<typename Scalar_, typename SizeType_>
    class AbstractVector
    // :
    //     public DistributedVector<Scalar_, SizeType_>,
    //     public Normed<Scalar_>,
    //     public Transformable<Scalar_>,
    //     public Reducible<Scalar_>,
    //     public Constructible<Scalar_, SizeType_, 1>,
    //     public ElementWiseOperand<Scalar_>,
    //     public ElementWiseOperand<AbstractVector<Scalar_, SizeType_>>,
    //     public Comparable<AbstractVector<Scalar_, SizeType_>>,
    //     public BLAS1Tensor<AbstractVector<Scalar_, SizeType_>>
        {
    public:
        using Scalar   = Scalar_;
        using SizeType = SizeType_;
        virtual ~AbstractVector() {}
    };

    template<class ConcreteType, int Order = Traits<ConcreteType>::Order>
    class Wrapper {};


    template<class Vector>
    class Wrapper<Vector, 1> : public AbstractVector<
                                        typename Traits<Vector>::Scalar,
                                        typename Traits<Vector>::SizeType
                                        > {
    public:

        template<class... Args>
        Wrapper(Args &&...args)
        : impl_( utopia::make_unique<Vector>(std::forward<Args>(args)...))
        {}

        template<class... Args>
        void construct(Args &&...args)
        {
            impl_ = utopia::make_unique<Vector>(std::forward<Args>(args)...);
        }

        std::unique_ptr<Vector> impl_;
    };
}

#endif //UTOPIA_ABSTRACT_VECTOR_HPP



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
    :
        public DistributedVector<Scalar_, SizeType_>
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

    template<class ConcreteType, int Order>
    class Traits<Wrapper<ConcreteType, Order> > : public Traits<ConcreteType> {};

    template<class Vector>
    class Wrapper<Vector, 1> : public AbstractVector<
                                        typename Traits<Vector>::Scalar,
                                        typename Traits<Vector>::SizeType
                                        > {
    public:
        using Scalar =  typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        template<class... Args>
        Wrapper(Args &&...args)
        : impl_( utopia::make_unique<Vector>(std::forward<Args>(args)...))
        {}

        template<class... Args>
        void construct(Args &&...args)
        {
            impl_ = utopia::make_unique<Vector>(std::forward<Args>(args)...);
        }

        inline SizeType local_size() const override
        {
            return impl_->local_size();
        }

        inline void c_set(const SizeType &i, const Scalar &value) override
        {
            impl_->c_set(i, value);
        }

        inline void c_add(const SizeType &i, const Scalar &value) override
        {
            impl_->c_add(i, value);
        }

        inline Range range() const override
        {
            return impl_->range();
        }

        //locks
        inline void read_lock()                override
        {
            impl_->read_lock();
        }

        inline void write_lock(WriteMode mode) override
        {
            impl_->write_lock(mode);
        }

        inline void read_unlock() override
        {
            impl_->read_unlock();
        }

        inline void write_unlock(WriteMode mode) override
        {
            impl_->write_unlock(mode);
        }

        inline void read_and_write_lock(WriteMode mode) override
        {
            impl_->read_and_write_lock(mode);
        }

        inline void read_and_write_unlock(WriteMode mode) override
        {
            impl_->read_and_write_unlock(mode);
        }

        //basic mutators
        inline void set(const SizeType &i, const Scalar &value) override
        {
            impl_->set(i, value);
        }

        inline void add(const SizeType &i, const Scalar &value) override
        {
            impl_->add(i, value);
        }

        inline Scalar get(const SizeType &i) const override
        {
            return impl_->get(i);
        }

        //print function
        inline void describe() const override
        {
            impl_->describe();
        }

        //utility functions
        inline bool empty() const override
        {
            return impl_->empty();
        }

        inline void clear() override
        {
            impl_->clear();
        }

        inline void set(const Scalar &val) override
        {
            impl_->set(val);
        }

        inline SizeType size() const override
        {
            return impl_->size();
        }

        inline Communicator &comm() override
        {
            return impl_->comm();
        }

        inline const Communicator &comm() const override
        {
            return impl_->comm();
        }

    private:
        std::unique_ptr<Vector> impl_;
    };
}

#endif //UTOPIA_ABSTRACT_VECTOR_HPP



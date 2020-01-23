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
    template<typename Scalar_, typename SizeType_>
    class Traits< AbstractVector<Scalar_, SizeType_> > {
    public:
        using Scalar = Scalar_;
        using SizeType = SizeType_;
    };

    //parallel types, collective operations
    template<typename Scalar_, typename SizeType_>
    class AbstractVector
    :
        public DistributedVector<Scalar_, SizeType_>,
        public Normed<Scalar_>,
        public Transformable<Scalar_>,
        public Reducible<Scalar_>,
        public Constructible<Scalar_, SizeType_, 1>,
        public ElementWiseOperand<Scalar_>,
        public ElementWiseOperand<AbstractVector<Scalar_, SizeType_>>,
        public Comparable<AbstractVector<Scalar_, SizeType_>>,
        public BLAS1Tensor<AbstractVector<Scalar_, SizeType_>>
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
        using Scalar         = typename Traits<Vector>::Scalar;
        using SizeType       = typename Traits<Vector>::SizeType;
        using AbstractVector = AbstractVector<
                                        typename Traits<Vector>::Scalar,
                                        typename Traits<Vector>::SizeType
                                        >;

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
        inline void read_lock() override
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
            std::cout << Traits<Vector>::backend_info().get_name() << " (vector)" << std::endl;
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

        inline Scalar norm_infty() const override
        {
            return impl_->norm_infty();
        }

        inline Scalar norm1() const override
        {
            return impl_->norm1();
        }

        inline Scalar norm2() const override
        {
            return impl_->norm2();
        }

        inline Scalar reduce(const Plus &op) const override
        {
            return impl_->reduce(op);
        }

        inline Scalar reduce(const Min &op)  const override
        {
            return impl_->reduce(op);
        }

        inline Scalar reduce(const Max &op)  const override
        {
            return impl_->reduce(op);
        }

        inline void e_mul(const Scalar &other) override
        {
            impl_->e_mul(other);
        }

        inline void e_div(const Scalar &other) override
        {
            impl_->e_div(other);
        }

        inline void e_min(const Scalar &other) override
        {
            impl_->e_min(other);
        }

        inline void e_max(const Scalar &other) override
        {
            impl_->e_max(other);
        }

        ////////////

        inline void e_mul(const AbstractVector &other) override
        {
            auto &other_w = static_cast<const Wrapper &>(other);
            assert(other_w.impl_);
            impl_->e_mul(*other_w.impl_);
        }

        inline void e_div(const AbstractVector &other) override
        {
            auto &other_w = static_cast<const Wrapper &>(other);
            assert(other_w.impl_);
            impl_->e_div(*other_w.impl_);
        }

        inline void e_min(const AbstractVector &other) override
        {
            auto &other_w = static_cast<const Wrapper &>(other);
            assert(other_w.impl_);
            impl_->e_min(*other_w.impl_);
        }

        inline void e_max(const AbstractVector &other) override
        {
            auto &other_w = static_cast<const Wrapper &>(other);
            assert(other_w.impl_);
            impl_->e_max(*other_w.impl_);
        }

        ///////////////////

        inline void swap(AbstractVector &x) override
        {
            auto &x_w = static_cast<const Wrapper &>(x);
            impl_->swap(*x_w.impl_);
        }

        ///<Scalar>SCAL - x = a*x
        inline void scale(const Scalar &a) override
        {
            impl_->scale(a);
        }

        ///<Scalar>COPY - copy x into y (this)
        inline void copy(const AbstractVector &x) override
        {
            auto &x_w = static_cast<const Wrapper &>(x);
            impl_->copy(*x_w.impl_);
        }

        ///<Scalar>AXPY - y = a*x + y
        inline void axpy(const Scalar &a, const AbstractVector &x) override
        {
            auto &x_w = static_cast<const Wrapper &>(x);
            impl_->axpy(a, *x_w.impl_);
        }

        ///<Scalar>DOT - dot product
        inline Scalar dot(const AbstractVector &x) const override
        {
            auto &x_w = static_cast<const Wrapper &>(x);
            return impl_->dot(*x_w.impl_);
        }

        inline bool equals(const AbstractVector &other, const Scalar &tol = 0.0) const override
        {
            auto &other_w = static_cast<const Wrapper &>(other);
            return impl_->equals(*other_w.impl_, tol);
        }

        //////////////////////

        inline void transform(const Sqrt &op) override
        {
            return impl_->transform(op);
        }

        inline void transform(const Pow2 &op) override
        {
            return impl_->transform(op);
        }

        inline void transform(const Log &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Exp &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Cos &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Sin &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Abs &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Minus &op) override
        {
            return impl_->transform(op);
        }

        inline void transform(const Pow &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Reciprocal<Scalar> &op) override
        {
            return impl_->transform(op);
        }

        inline void zeros(const SizeType &s) override { impl_->zeros(s); }
        inline void values(const SizeType &s, const Scalar &val) override { impl_->values(s, val); }

        inline void local_zeros(const SizeType &s) override { impl_->local_zeros(s); }
        inline void local_values(const SizeType &s, const Scalar &val) override { impl_->local_values(s, val); }

        //comodity
        inline void zeros(const Size &s) override { impl_->zeros(s); }
        inline void values(const Size &s, const Scalar &val) override { impl_->values(s, val); }

        inline void local_zeros(const Size &s) override { impl_->local_values(s, 0.0); }
        inline void local_values(const Size &s, const Scalar &val) override { impl_->local_values(s, val); }


    private:
        std::unique_ptr<Vector> impl_;
    };
}

#endif //UTOPIA_ABSTRACT_VECTOR_HPP



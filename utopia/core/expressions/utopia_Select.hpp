#ifndef UTOPIA_SELECT_HPP
#define UTOPIA_SELECT_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Utils.hpp"

namespace utopia {
    template<class Expr, int Order>
    class Select;

    template<class Expr_>
    class Select<Expr_, 1> : public Expression< Select<Expr_, 1>  > {
    public:
        typedef Expr_ Expr;
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        //FIXME use Traits instead
        static const int Order = 1;
        static_assert(Expr::Order == Order, "must be same order of the tensor");

        inline explicit Select(const Expr &expr, const std::vector<SizeType> &index)
        : expr_(expr), index_ptr_(utopia::make_ref(index))
        {}

        inline explicit Select(const Expr &expr, std::vector<SizeType> &&index)
        : expr_(expr), index_ptr_(std::make_shared(std::move(index)))
        {}

        inline const std::vector<SizeType> &index() const
        {
            return *index_ptr_;
        }

        const Expr &expr() const
        {
            return expr_;
        }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
        std::shared_ptr<const std::vector<SizeType> > index_ptr_;
    };


    template<class Expr_>
    class Select<Expr_, 2> : public Expression< Select<Expr_, 2> > {
    public:
        typedef Expr_ Expr;
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        //FIXME use Traits instead
        static const int Order = 2;
        static_assert(Expr::Order == Order, "must be same order of the tensor");


        inline explicit Select(const Expr &expr, const std::vector<SizeType> &row_index, const std::vector<SizeType> &col_index)
        : expr_(expr),
          row_index_ptr_(utopia::make_ref(row_index)),
          col_index_ptr_(utopia::make_ref(col_index))
        {}

        inline explicit Select(const Expr &expr, std::vector<SizeType> &&row_index, std::vector<SizeType> &&col_index)
        : expr_(expr),
          row_index_ptr_(std::make_shared(std::move(row_index))),
          col_index_ptr_(std::make_shared(std::move(col_index)))
        {}

        inline const std::vector<SizeType> &row_index() const
        {
            return *row_index_ptr_;
        }

        inline const std::vector<SizeType> &col_index() const
        {
            return *col_index_ptr_;
        }

        const Expr &expr() const
        {
            return expr_;
        }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
        std::shared_ptr<const std::vector<SizeType> > row_index_ptr_;
        std::shared_ptr<const std::vector<SizeType> > col_index_ptr_;
    };

    template<class Expr, int Order>
    class Traits< Select<Expr, Order> > : public Traits<Expr> {};


    template<class Derived>
    class Selectable {};

    template<class Derived>
    class Selectable<Tensor<Derived, 1>> {
    public:
        using TensorT = utopia::Tensor<Derived, 1>;
        using That    = utopia::Selectable<Tensor<Derived, 1>>;

        typedef typename utopia::Traits<Derived>::SizeType SizeType;

        //lazy evaluation
        inline friend Select<TensorT, 1> select(const That &that, const std::vector<SizeType> &index)
        {
            return Select<TensorT, 1>(that.derived(), index);
        }

        inline friend Select<TensorT, 1> select(const That &that, std::vector<SizeType> &&index)
        {
            return Select<TensorT, 1>(that.derived(), std::move(index));
        }

        //direct evaluation
        virtual void select(const std::vector<SizeType> &index, Derived &result) const = 0;

    private:
        CONST_DERIVED_CRT(TensorT);
    };

    template<class Derived>
    class Selectable<Tensor<Derived, 2>> {
    public:
        using TensorT = utopia::Tensor<Derived, 2>;
        using That = utopia::Selectable<Tensor<Derived, 2>>;

        typedef typename utopia::Traits<Derived>::SizeType SizeType;

        //lazy evaluation
        inline friend Select<TensorT, 2> select(const That &that, const std::vector<SizeType> &row_index, const std::vector<SizeType> &col_index = std::vector<SizeType>())
        {
            return Select<TensorT, 2>(that.derived(), row_index, col_index);
        }

        inline friend Select<TensorT, 2> select(const That &that, std::vector<SizeType> &&row_index, std::vector<SizeType> &&col_index = std::vector<SizeType>())
        {
            return Select<TensorT, 2>(that.derived(), std::move(row_index, col_index));
        }

        /// if col_index is empty select all columns
        //direct evaluation
        virtual void select(
            const std::vector<SizeType> &row_index, 
            const std::vector<SizeType> &col_index, 
            Derived &result) const = 0;

    private:
        CONST_DERIVED_CRT(Derived);
    };

}

#endif //UTOPIA_SELECT_HPP

#ifndef UTOPIA_SELECT_HPP
#define UTOPIA_SELECT_HPP 

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Utils.hpp"

namespace utopia {
	template<class Expr, typename SizeType, int Order>
	class Select;

	template<class Expr_, typename SizeType>
	class Select<Expr_, SizeType, 1> : public Expression< Select<Expr_, SizeType, 1> > {
	public:

		typedef Expr_ Expr;
		typedef typename Expr::Scalar Scalar;
		
		enum {
		    Order = Expr::Order
		};


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


	template<class Expr_, typename SizeType>
	class Select<Expr_, SizeType, 2> : public Expression< Select<Expr_, SizeType, 2> > {
	public:
		typedef Expr_ Expr;
		typedef typename Expr::Scalar Scalar;

		enum {
		    Order = Expr::Order
		};

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

	template<class Expr, typename SizeType, int Order>
	class Traits< Select<Expr, SizeType, Order> > : public Traits<Expr> {};


	template<class Implementation, class Derived, int Order>
	class Selectable {};

	template<class Implementation, class Derived>
	class Selectable<Implementation, Derived, 1> {
	public:
		typedef typename utopia::Traits<Implementation>::SizeType SizeType;

	    inline Select<Derived, SizeType, 1> select(const std::vector<SizeType> &index) const
	    {
	        return Select<Derived, SizeType, 1>(derived(), index);
	    }

	    inline Select<Derived, SizeType, 1> select(std::vector<SizeType> &&index) const
	    {
	        return Select<Derived, SizeType, 1>(derived(), std::move(index));
	    }

	private:
	    CONST_DERIVED_CRT(Derived);
	};


	template<class Implementation, class Derived>
	class Selectable<Implementation, Derived, 2> {
	public:
		typedef typename utopia::Traits<Implementation>::SizeType SizeType;

	    inline Select<Derived, SizeType, 2> select(const std::vector<SizeType> &row_index, const std::vector<SizeType> &col_index = std::vector<SizeType>()) const
	    {
	        return Select<Derived, SizeType, 2>(derived(), row_index, col_index);
	    }

	    inline Select<Derived, SizeType, 2> select(std::vector<SizeType> &&row_index, std::vector<SizeType> &&col_index = std::vector<SizeType>()) const
	    {
	        return Select<Derived, SizeType, 2>(derived(), std::move(row_index, col_index));
	    }

	private:
	    CONST_DERIVED_CRT(Derived);
	};

}

#endif //UTOPIA_SELECT_HPP

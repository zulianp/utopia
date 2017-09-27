#ifndef UTOPIA_SELECT_HPP
#define UTOPIA_SELECT_HPP 

#include "utopia_Utils.hpp"

namespace utopia {
	template<class Expr, typename SizeType, int Order>
	class Select;

	template<class Expr, typename SizeType>
	class Select<Expr, SizeType, 1> : public Expression< Select<Expr, SizeType, 1> > {
	public:
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

	template<class Expr, typename SizeType, int Order>
	class Traits< Select<Expr, SizeType, Order> > : public Traits<Expr> {};


	// template<class Derived, int Order>
	// class Selectable {};

	// template<class Derived>
	// class Selectable<Derived, 2> {
	// public:
	//     Select<Derived, 2> select(const std::vector<SizeType> &row_index, const std::vector<SizeType> &col_index)
	//     {
	//         return Select<Derived, 2>(derived(), row_index, col_index);
	//     }

	// private:
	//     CONST_DERIVED_CRT(Derived);
	// };

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

}

#endif //UTOPIA_SELECT_HPP

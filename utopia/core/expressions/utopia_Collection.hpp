// #ifndef UTOPIA_ITERABLE_HPP
// #define UTOPIA_ITERABLE_HPP 

// #include <iterator>
// #include "utopia_ForwardDeclarations.hpp"
// #include "utopia_Assign.hpp"

// namespace utopia {
// 	template<class Expr>
// 	class Collection : public Expression< Collection<Expr> > {
// 	public:
// 			typedef utopia::Traits<Expr> Traits;
// 			typedef typename TypeAndFill<Traits, Expr>::Type Result;
// 			static const int Order = Traits::Order;

// 		Collection(std::vector<Expr> &container)
// 		: container_(utopia::make_ref(container))
// 		{}

// 		Collection()
// 		: container_(std::make_shared< std::vector<Expr> >())
// 		{}

// 		inline operator std::vector< Wrapper<Result, Order> >() const
// 		{
// 			std::vector< Wrapper<Result, Order> > result(get().size());

// 			Evaluator<typename Traits::Vector, Traits::Backend> eval;

// 			for(std::size_t i = 0; i < get().size(); ++i) {
// 				eval.eval( Construct< Wrapper<Result, Order>, Expr >( result[i], get()[i] ) );
// 			}

// 			return result;
// 		}

// 		std::vector<Expr> &get()
// 		{
// 			return *container_;
// 		}

// 		const std::vector<Expr> &get() const
// 		{
// 			return *container_;
// 		}


// 	private:
// 		std::shared_ptr<std::vector<Expr> > container_;
// 	};

// 	template<class Expr>       
// 	Collection<Expr> wrap(std::vector<Expr> &container)
// 	{
// 		return Collection<Expr>(container);
// 	}

// 	template<class Left, class Right>
// 	Collection<Multiply<Left, Right> > operator*(const Left &left, const Collection< Wrapper<Right, Order> > &right)
// 	{
// 		Collection<Multiply<Left, Right> > ret;
// 		ret.get().resize(right.get().size());

// 		for(std::size_t i = 0; i < right.get().size(); ++i) {
// 			ret.get()[i] = left * right.get()[i];
// 		}

// 		return ret;
// 	}

// 	template<class Expr>
// 	class Traits< Collection<Expr> > : public Traits<Expr> {};
// }

// #endif //UTOPIA_ITERABLE_HPP

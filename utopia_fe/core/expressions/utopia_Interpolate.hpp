#ifndef UTOPIA_FE_INTERPOLATE_HPP
#define UTOPIA_FE_INTERPOLATE_HPP 

#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_Any.hpp"

#include <utility>

namespace utopia {

	template<class Coefficient, class Fun>
	class Interpolate : public Expression< Interpolate<Coefficient, Fun> > {
	public:

		typedef typename Traits<Coefficient>::Scalar Scalar;

		enum {
			Order = Traits<Coefficient>::Order
		};

		template<class Coeff>
		Interpolate(Coeff &&coeff, Fun &fun)
		: coeff_(std::forward<Coeff>(coeff)), fun_(fun)
		{}

		inline Fun &fun()
		{
			return fun_;
		}

		inline const Fun &fun() const
		{
			return fun_;
		}

		inline Coefficient &coefficient()
		{
			return coeff_;
		}

		inline const Coefficient &coefficient() const
		{
			return coeff_;
		}

		inline std::string getClass() const override {
			return "Interpolate";
		}

	private:
		UTOPIA_STORE(Coefficient) coeff_;
		UTOPIA_STORE(Fun) fun_;	
	};

	template<>
	class Interpolate<Any, Any> {};

	template<class Coefficient, class Fun>
	class Traits< Interpolate<Coefficient, Fun> > : public Traits<Coefficient> {
	public:
	
	};


	template<class Coefficient, class Fun>
	inline Interpolate<Coefficient, Fun> interpolate(Coefficient &coeff, Fun &fun)
	{
		return Interpolate<Coefficient, Fun>(coeff, fun);
	}


	template<class Coefficient, class Fun, class Tensor>
	inline Interpolate<Coefficient, Fun> interpolate(Coefficient &coeff, Fun &fun, Tensor &&tensor)
	{
		return Interpolate<Coefficient, Fun>(coeff, fun, tensor);
	}

	//FIXME: Deprecated remove once refactor is finished
	template<class T, int Order, class Fun, class Tensor>
	inline Interpolate<ConstantCoefficient<T, Order>, Fun> interpolate(const ConstantCoefficient<T, Order> &coeff, Fun &fun, Tensor &&tensor)
	{
		return Interpolate<ConstantCoefficient<T, Order>, Fun>(coeff, fun, tensor);
	}
	
}

#endif //UTOPIA_FE_INTERPOLATE_HPP

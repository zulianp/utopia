#ifndef UTOPIA_FE_INTERPOLATE_HPP
#define UTOPIA_FE_INTERPOLATE_HPP 

#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"


namespace utopia {

	template<class Coefficient, class Fun>
	class Interpolate : public Expression< Interpolate<Coefficient, Fun> > {
	public:

		typedef typename Traits<Coefficient>::Scalar Scalar;

		enum {
			Order = Traits<Coefficient>::Order
		};

		Interpolate(const Coefficient &coeff, Fun &fun)
		: coeff_(coeff), fun_(fun)
		{}


	private:
		Coefficient coeff_;
		UTOPIA_STORE(Fun) fun_;	
	};

	template<class Coefficient, class Fun>
	class Traits< Interpolate<Coefficient, Fun> > : public Traits<Coefficient> {
	public:
	
	};

		

	template<class Coefficient, class Fun>
	inline Interpolate<Coefficient, Fun> interpolate(const Coefficient &coeff, Fun &fun)
	{
		return Interpolate<Coefficient, Fun>(coeff, fun);
	}


	template<class Coefficient, class Fun, class Tensor>
	inline Interpolate<Coefficient, Fun> interpolate(const Coefficient &coeff, Fun &fun, Tensor &&tensor)
	{
		return Interpolate<Coefficient, Fun>(coeff, fun, tensor);
	}
	
}

#endif //UTOPIA_FE_INTERPOLATE_HPP

#ifndef UTOPIA_FE_INTERPOLATE_HPP
#define UTOPIA_FE_INTERPOLATE_HPP 

#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"


namespace utopia {

	template<class Coefficient, class FESpace>
	class Interpolate : public Expression< Interpolate<Coefficient, FESpace> > {
	public:

		typedef typename Traits<Coefficient>::Scalar Scalar;

		enum {
			Order = Traits<Coefficient>::Order
		};

		Interpolate(const Coefficient &coeff, FESpace &fe)
		: coeff_(coeff), fe_(fe)
		{}


	private:
		Coefficient coeff_;
		UTOPIA_STORE(FESpace) fe_;	
	};

	template<class Coefficient, class FESpace>
	class Traits< Interpolate<Coefficient, FESpace> > : public Traits<Coefficient> {
	public:
	
	};

		

	template<class Coefficient, class FESpace>
	inline Interpolate<Coefficient, FESpace> interpolate(const Coefficient &coeff, FESpace &fe)
	{
		return Interpolate<Coefficient, FESpace>(coeff, fe);
	}


	template<class Coefficient, class FESpace, class Tensor>
	inline Interpolate<Coefficient, FESpace> interpolate(const Coefficient &coeff, FESpace &fe, Tensor &&tensor)
	{
		return Interpolate<Coefficient, FESpace>(coeff, fe, tensor);
	}
	
}

#endif //UTOPIA_FE_INTERPOLATE_HPP

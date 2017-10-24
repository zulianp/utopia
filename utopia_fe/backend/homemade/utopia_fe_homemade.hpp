#ifndef UTOPIA_FE_HOMEMADE_HPP
#define UTOPIA_FE_HOMEMADE_HPP 


#include "utopia_Traits.hpp"
#include "utopia_Base.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_intersector.hpp"

namespace utopia {

	static const int EMPTY_BACKEND_FLAG = -3;

	class Empty {};

	template<>
	class Traits<Empty> {
	public:
		static const int Backend = EMPTY_BACKEND_FLAG;
		static const int Order = 1;
		typedef double Scalar;

	};


	class EmptyFunctionSpace : public FunctionSpace<EmptyFunctionSpace> {
	public:
	};

	template<>
	class FormTraits<EmptyFunctionSpace> : public Traits<Empty> {
	public:
		typedef utopia::Empty Implementation;
	};

	template<typename T>
	class Elem {
	public:
		//geom shapes
		static const int TRIANGLE = 0;

		T * points() { return nullptr; }
		int n_nodes() { return 0; }
		int type() { return TRIANGLE; }
	};

	template<typename T>
	class Quadrature {
	public:
		T * points() { return nullptr; }
		int n_points() { return 0; }
	};

	template<typename T, int FEType, int FEOrder>
	class FE {
	public:
		//fe types
		static const int LAGRANGE = 0;

		void init(Elem<T> &elem, Quadrature<T> &quad)
		{
			switch(elem.type()) {
				case Elem<T>::TRIANGLE:
				{
					isect.triangle_make_fe_object(elem.points(), quad.n_points(), quad.points(), &fe_);
					break;
				}
				default: 
				{
					assert(false);
					break;
				}
			}
		}

	private:
		Intersector isect;
		FEObject fe_;
	};
}

#endif //UTOPIA_FE_HOMEMADE_HPP

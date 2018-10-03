#ifndef UTOPIA_KOKKOS_OP_HPP
#define UTOPIA_KOKKOS_OP_HPP
#include <Kokkos_Core.hpp>
#include "utopia_Operators.hpp"
#include <iostream>


namespace utopia {
	template<class Scalar, class Op> class KokkosOp;

		template<class Scalar>
		class KokkosOp<Scalar, Max> {
		public:
		KOKKOS_FUNCTION	static double apply(const Scalar &a, const  Scalar &b) {
	        	if(a > b){
	        		return a;
	        	}
	        	else
	        	{
	        		return b;
	        	}
	        }
		};

		template<class Scalar>
		class KokkosOp<Scalar, Min> {
		public:
		KOKKOS_FUNCTION	static double apply(const Scalar &a, const  Scalar &b) {
	        	if(a > b){
	        		return b;
	        	}
	        	else
	        	{
	        		return a;
	        	}
	        }
		};


}
#endif 



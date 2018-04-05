#ifndef UTOPIA_NONLINEAR_SOLVERS_INTERFACES_HPP
#define UTOPIA_NONLINEAR_SOLVERS_INTERFACES_HPP 

#include "utopia_ForwardDeclarations.hpp"

namespace utopia 
{

  	template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
  	class NonLinearGaussSeidel 
  	{
  		public: 
  			NonLinearGaussSeidel()
  			{
  				static_assert(Backend < HOMEMADE, "NonLinearGaussSeidel not implemented for this backend");
  			}
  	};


  	template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
  	class NonLinearConjugateGradient 
  	{
  		public: 
  			NonLinearConjugateGradient()
  			{
  				static_assert(Backend < HOMEMADE, "NonLinearConjugateGradient not implemented for this backend");
  			}
  	};

  	template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
  	class NonLinearGMRES
  	{
  		public: 
  			NonLinearGMRES()
  			{
  				static_assert(Backend < HOMEMADE, "NonLinearGMRES not implemented for this backend");
  			}
  	};



  	template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
  	class NonLinearAnderson
  	{
  		public: 
  			NonLinearAnderson()
  			{
  				static_assert(Backend < HOMEMADE, "NonLinearAnderson not implemented for this backend");
  			}
  	};


	template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
	class NonLinearRichardson
	{
  		public: 
	  		NonLinearRichardson()
	  		{
	  			static_assert(Backend < HOMEMADE, "NonLinearAnderson not implemented for this backend");
	  		}
  	};


  	// line-search version
  	// template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
  	// class Newton
  	// {
  	// 	public: 
  	// 		Newton()
  	// 		{
  	// 			static_assert(Backend < HOMEMADE, "NonLinearAnderson not implemented for this backend");
  	// 		}
  	// };


   //  template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
  	// class TrustRegion
  	// {
  	// 	public: 
	  // 		TrustRegion()
	  // 		{
	  // 			static_assert(Backend < HOMEMADE, "NonLinearAnderson not implemented for this backend");
	  // 		}
  	// };


}

#endif //UTOPIA_NONLINEAR_SOLVERS_INTERFACES_HPP

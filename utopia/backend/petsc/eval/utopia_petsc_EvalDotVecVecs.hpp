#ifndef UTOPIA_EVAL_PETSC_DOT_VEC_VECS_HPP
#define UTOPIA_EVAL_PETSC_DOT_VEC_VECS_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia 
{

    template<class Vector>
    class EvalDots <Vector, PETSC>
    {
    	public:
	        static void apply(const Wrapper<Vector, 1> &v1, const std::vector<std::shared_ptr<Wrapper<Vector, 1> > > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results) 
	        { 
	        	typename utopia::Traits<Vector>::SizeType n =  vectors.size(); 

	        	// array of the dot products should not be allocated
		    	if(n!=results.size())
		    		results.resize(n); 

		    	std::vector<Vec> vecs(n);

		    	for(auto i=0; i < vectors.size(); i++)
		    		vecs[i]=(raw_type(*vectors[i])); 

	        	 VecMDot(v1.implementation().implementation(), n, vecs.data(), results.data()); 
	        }
    };

}


#endif //UTOPIA_EVAL_PETSC_DOT_VEC_VECS_HPP

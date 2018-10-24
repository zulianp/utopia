#ifndef UTOPIA_EVAL_DOT_VEC_VECS_HPP
#define UTOPIA_EVAL_DOT_VEC_VECS_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia 
{
    template<class Vector, int Backend = Traits<Vector>::Backend>
    class EvalDots 
    {
    	public:
	        static void apply(const Wrapper<Vector, 1> &v1, const std::vector<std::shared_ptr<Wrapper<Vector, 1> > > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results) 
	        { 
	        	typename utopia::Traits<Vector>::SizeType n =  vectors.size(); 

		    	if(n!=results.size())
		    		results.resize(n); 

		    	for(auto i = 0; i < n; i++)
		    		results[i] = dot(v1, *(vectors[i])); 
	        }

	        static void apply(const Wrapper<Vector, 1> &v1, const std::vector<Wrapper<Vector, 1> > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results) 
	        { 
	        	typename utopia::Traits<Vector>::SizeType n =  vectors.size(); 

		    	if(n!=results.size())
		    		results.resize(n); 

		    	for(auto i = 0; i < n; i++)
		    		results[i] = dot(v1, (vectors[i])); 
	        }	        
    };


	template< class Vector>
	void dots(const Wrapper<Vector, 1> &v1, const std::vector<std::shared_ptr<Wrapper<Vector, 1> > > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results) 
	{ 
		EvalDots<Vector>::apply(v1, vectors, results); 
	}


	template< class Vector>
	void dots(const Wrapper<Vector, 1> &v1, const std::vector<Wrapper<Vector, 1> > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results) 
	{ 
		EvalDots<Vector>::apply(v1, vectors, results); 
	}


}


#endif //UTOPIA_EVAL_DOT_VEC_VECS_HPP

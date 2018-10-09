#ifndef UTOPIA_MAT_GET_COL_HPP
#define UTOPIA_MAT_GET_COL_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia 
{
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class EvalGetCol 
    {
    	public:
	        static void apply(const Wrapper<Matrix, 2> &M, Wrapper<Vector, 1> &v, typename utopia::Traits<Vector>::SizeType col_id)
	        { 
		     	using VectorT  = utopia::Wrapper<Vector, 1>;
		        using MatrixT  = utopia::Wrapper<Matrix, 2>;

		        using SizeType = UTOPIA_SIZE_TYPE(VectorT);        

		        SizeType local_rows = local_size(M).get(0); 
		        SizeType global_cols = size(M).get(1); 

		        if(global_cols< col_id)
		            std::cerr<<"mat_get_col: Requested column id does not exist. \n"; 

		        if(local_rows!= local_size(v).get(0) || empty(v))
		            v = local_zeros(local_rows); 

		        {
		            Read<MatrixT> r_m(M);
		            Write<VectorT> w_v(v);

		            auto r = row_range(M);

		            for(auto i = r.begin(); i != r.end(); ++i) 
		            {
		                v.set(i, M.get(i, col_id));
		            }
		        }
	        }
    };


	template<class Matrix, class Vector>
	void mat_get_col(const Wrapper<Matrix, 2> &M, Wrapper<Vector, 1> &v, typename utopia::Traits<Vector>::SizeType col_id)
	{ 
		EvalGetCol<Matrix, Vector>::apply(M, v, col_id); 
	}



}


#endif //UTOPIA_MAT_GET_COL_HPP

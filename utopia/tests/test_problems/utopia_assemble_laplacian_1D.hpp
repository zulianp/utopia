#ifndef UTOPIA_ASSEMBLE_LAPLACIAN_1D
#define UTOPIA_ASSEMBLE_LAPLACIAN_1D

#include "utopia_Base.hpp"
#include "utopia_Range.hpp"
#include "utopia_Writable.hpp"

namespace utopia {

	template<class Matrix>
	void assemble_laplacian_1D(const utopia::SizeType n, Matrix &m)
	{

	    // n x n matrix with maximum 3 entries x row        
		{
			Write<Matrix> w(m);
			Range r = row_range(m);

	        //You can use set instead of add. [Warning] Petsc does not allow to mix add and set.
			for(SizeType i = r.begin(); i != r.end(); ++i) {
				if(i > 0) {    
					m.add(i, i - 1, -1.0);    
				}

				if(i < n-1) {
					m.add(i, i + 1, -1.0);
				}

				m.add(i, i, 2.0);
			}
		}
	}
}

#endif //UTOPIA_ASSEMBLE_LAPLACIAN_1D

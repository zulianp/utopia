#ifndef UTOPIA_ASSEMBLE_LAPLACIAN_1D
#define UTOPIA_ASSEMBLE_LAPLACIAN_1D

#include "utopia_Base.hpp"
#include "utopia_Range.hpp"
#include "utopia_Writable.hpp"

namespace utopia {

	template<class Matrix>
	void assemble_laplacian_1D(Matrix &m, const bool bc = false)
	{
	    // n x n matrix with maximum 3 entries x row        
		Write<Matrix> w(m);
		Range r = row_range(m);
		auto n = size(m).get(0);
		
		for(SizeType i = r.begin(); i != r.end(); ++i) {
			if(bc) {
				if(i == 0) {
					m.set(0, 0, 1.);
					continue;
				}

				if(i == n-1) {
					m.set(n-1, n-1, 1.);
					continue;
				}
			}

			if(i > 0) {    
				m.set(i, i - 1, -1.0);    
			}

			if(i < n-1) {
				m.set(i, i + 1, -1.0);
			}

			if(i == 0 || i == n-1) {
				m.set(i, i, 1.);
			} else {
				m.set(i, i, 2.0);
			}
		}
	}
}

#endif //UTOPIA_ASSEMBLE_LAPLACIAN_1D

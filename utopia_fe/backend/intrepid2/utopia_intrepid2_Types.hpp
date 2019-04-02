#ifndef UTOPIA_INTREPID2_TYPES_HPP
#define UTOPIA_INTREPID2_TYPES_HPP

#include <iostream>

//#include "utopia_fe_kokkos_fix.hpp"

#include "utopia_intrepid2_FEForwardDeclarations.hpp"
#include "utopia.hpp"
#include "utopia_fe_base.hpp"

/*#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"*/

namespace utopia {

	inline static void add_matrix(const std::vector<libMesh::dof_id_type> &row_dofs,
								  const std::vector<libMesh::dof_id_type> &col_dofs,
								  TSMatrixd &mat) //USparseMatrix
	{
		int size = row_dofs.size()*col_dofs.size();
		std::vector<double> zero_vector (size);
		for(int i=0;i<size;i++)
			zero_vector[i]=0.;

		mat.add_matrix<libMesh::dof_id_type>(row_dofs, col_dofs, zero_vector);
		//mat.add_matrix<libMesh::dof_id_type>(row_dofs, col_dofs, 0.);
	}


}

#endif //UTOPIA_INTREPID2_TYPES_HPP
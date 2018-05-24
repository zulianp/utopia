#ifndef UTOPIA_LOCAL_2_GLOBAL_HPP
#define UTOPIA_LOCAL_2_GLOBAL_HPP

#include "moonolith_communicator.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_sparse_matrix.hpp"

#include "libmesh/dense_matrix.h"

namespace utopia {

	class Local2Global {
	public:
		using LocalMatrix = libMesh::DenseMatrix<libMesh::Real>;

		Local2Global(const bool use_set_instead_of_add)
		: use_set_instead_of_add(use_set_instead_of_add)
		{}

		void apply(
			const std::vector<long> &trial,
			const std::vector<long> &test,
			const LocalMatrix &local_mat,
			moonolith::SparseMatrix<double> &global_mat) const;
		
		void add(
			const std::vector<long> &trial,
			const std::vector<long> &test,
			const LocalMatrix &local_mat,
			moonolith::SparseMatrix<double> &global_mat) const;

		void set_non_zero(
			const std::vector<long> &trial,
			const std::vector<long> &test,
			const LocalMatrix &local_mat,
			moonolith::SparseMatrix<double> &global_mat) const;

		void redistribute(
			moonolith::Communicator &comm,
			const long n_local_dofs_trial,
			const long n_local_dofs_test,
			moonolith::SparseMatrix<double> &global_mat);

		bool check_valid_matrix(const moonolith::SparseMatrix<double> &global_mat) const;

	private:
		bool use_set_instead_of_add;
	};

}

#endif //UTOPIA_LOCAL_2_GLOBAL_HPP

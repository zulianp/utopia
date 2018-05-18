#include "utopia_Local2Global.hpp"

#include <numeric>

namespace utopia {
	void Local2Global::apply(
		const std::vector<long> &trial,
		const std::vector<long> &test,
		const LocalMatrix &local_mat,
		moonolith::SparseMatrix<double> &global_mat) const
	{
		if(use_set_instead_of_add) {
			set_non_zero(trial, test, local_mat, global_mat);
		} else {
			add(trial, test, local_mat, global_mat);
		}
	}

	void Local2Global::add(
		const std::vector<long> &trial,
		const std::vector<long> &test,
		const LocalMatrix &local_mat,
		moonolith::SparseMatrix<double> &global_mat) const
	{
		for(std::size_t i = 0; i < test.size(); ++i) {
			const auto dof_I = test[i];

			for(std::size_t j = 0; j < trial.size(); ++j) {
				const auto dof_J = trial[j];

				global_mat.add(dof_I, dof_J, local_mat(i, j));
			}
		}
	}

	void Local2Global::set_non_zero(
		const std::vector<long> &trial,
		const std::vector<long> &test,
		const LocalMatrix &local_mat,
		moonolith::SparseMatrix<double> &global_mat) const
	{
		for(std::size_t i = 0; i < test.size(); ++i) {
			const auto dof_I = test[i];

			for(std::size_t j = 0; j < trial.size(); ++j) {
				const auto dof_J = trial[j];

				if(std::abs(local_mat(i, j)) != 0.) {
					global_mat.set(dof_I, dof_J, local_mat(i, j));

					//REMOVE ME
					// assert(global_mat.at(dof_I, dof_J) < 1.0001);
				}
			}
		}
	}

	void Local2Global::redistribute(
		moonolith::Communicator &comm,
		const long n_local_dofs_trial,
		const long n_local_dofs_test,
		moonolith::SparseMatrix<double> &global_mat)
	{
		std::vector<moonolith::Integer> range_master(comm.size() + 1, 0);
		std::vector<moonolith::Integer> range_slave(comm.size()  + 1, 0);

		range_master[comm.rank() + 1] += static_cast<unsigned int>(n_local_dofs_trial);
		range_slave [comm.rank() + 1] += static_cast<unsigned int>(n_local_dofs_test);

		comm.all_reduce(&range_master[0], range_master.size(), moonolith::MPISum());
		comm.all_reduce(&range_slave[0],  range_slave.size(),  moonolith::MPISum());

		std::partial_sum(range_master.begin(), range_master.end(), range_master.begin());
		std::partial_sum(range_slave.begin(),  range_slave.end(),  range_slave.begin());

		moonolith::Redistribute< moonolith::SparseMatrix<double> > redist(comm.get_mpi_comm());

		if(use_set_instead_of_add) {
			// assert(check_valid_matrix(global_mat));
			
			redist.apply(range_slave, global_mat, moonolith::Assign<double>());

			// assert(check_valid_matrix(global_mat));
		} else {
			redist.apply(range_slave, global_mat, moonolith::AddAssign<double>());
		}

		assert(range_slave.empty() == range_master.empty() || range_master.empty());
	}


	bool Local2Global::check_valid_matrix(const moonolith::SparseMatrix<double> &global_mat) const
	{
		auto vec = global_mat.sum_cols();

		for(auto v : vec) {
			if(v > 1.0001) {
				assert(false);
				return false;
			}
		}

		return true;
	}
}

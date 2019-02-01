#include "utopia_private_PetrovGalerkinAssembler.hpp"

namespace utopia {

	namespace private_ {

		PetrovGalerkinAssembler::PetrovGalerkinAssembler()
		{}

		bool PetrovGalerkinAssembler::assemble(
			const Elem &master,
			FEType master_type,
			const Elem &slave,
			FEType slave_type,
			std::function<void(std::vector<long> &, std::vector<long> &)> dof_fun
			)
		{
			for(auto &mat_i : elemmat_) {
				mat_i.zero();
			}

			if(assembler_->assemble(master,
				master_type,
				slave,
				slave_type,
				elemmat_)) {

				dof_fun(master_dofs_, slave_dofs_);

			++n_intersections_;

			for(std::size_t i = 0; i < elemmat_.size(); ++i) {
				auto &mat_i = elemmat_[i];
				if(mat_i.get_values().empty()) continue;

				auto partial_sum = std::accumulate(mat_i.get_values().begin(), mat_i.get_values().end(), libMesh::Real(0.0));

				assert(!std::isnan(partial_sum));
				local_element_matrices_sum_[i] += partial_sum;

				assert(static_cast<bool>(mat_buffer_[i]));

				switch(assembler_->type(i)) {

					case LocalAssembler::MASTER_X_SLAVE:
					{
						local2global_->apply(master_dofs_, slave_dofs_, elemmat_[i], *mat_buffer_[i]);
						break;
					}

					case LocalAssembler::SLAVE_X_SLAVE:
					{
						local2global_->apply(slave_dofs_, slave_dofs_, elemmat_[i], *mat_buffer_[i]);
						break;
					}

					case LocalAssembler::MASTER_X_MASTER:
					{
						local2global_->apply(master_dofs_, master_dofs_, elemmat_[i], *mat_buffer_[i]);
						break;
					}

					case LocalAssembler::SLAVE_X_MASTER:
					{
						local2global_->apply(slave_dofs_, master_dofs_, elemmat_[i], *mat_buffer_[i]);
						break;
					}

					default:
					{
						assert(false);
						break;
					}
				}
			}

			return true;
		} else {
			n_false_positives_++;
			return false;
		}
	}

	void PetrovGalerkinAssembler::initialize(
		const moonolith::Communicator &comm,
		const std::shared_ptr<LocalAssembler> &assembler,
		const std::shared_ptr<Local2Global> &local2global,
		const TransferOptions &opts,
		const SizeType from_n_dofs,
		const SizeType from_n_local_dofs,
		const SizeType to_n_dofs,
		const SizeType to_n_local_dofs)
	{
		comm_ = comm;
		assembler_ = assembler;
		local2global_ = local2global;
		opts_ = opts;

		from_n_dofs_ = from_n_dofs;
		from_n_local_dofs_ = from_n_local_dofs;
		to_n_dofs_ = to_n_dofs;
		to_n_local_dofs_ = to_n_local_dofs;

		n_intersections_ = 0;
		n_false_positives_ = 0;

		init_buffers();
	}

	void PetrovGalerkinAssembler::finalize(
		std::vector<std::shared_ptr<PetrovGalerkinAssembler::SparseMatrix>> &mats)
	{
		auto n_forms = assembler_->n_forms();
		mats.resize(n_forms);

		for(std::size_t i = 0; i < n_forms; ++i) {
			if(!mats[i]) {
				mats[i] = std::make_shared<SparseMatrix>();
			}

			SparseMatrix &mat = *mats[i];
			finalize_form(i, mat);
		}
	}

	bool PetrovGalerkinAssembler::check_n_forms(const int n_forms)
	{
		int n_forms_max = n_forms;
		comm_.all_reduce(&n_forms_max, 1, moonolith::MPIMax());
		bool ok = n_forms == n_forms_max;

		if(!ok) {
			std::cerr << comm_ << n_forms << " != " << n_forms_max << std::endl;
		}

		assert(ok);
		return ok;
	}


	void PetrovGalerkinAssembler::init_buffers() {
		assert(assembler_);

		auto n_forms = assembler_->n_forms();
		assert(n_forms > 0);

		assert(check_n_forms(n_forms));

		mat_buffer_.resize(n_forms);
		elemmat_.resize(n_forms);
		local_element_matrices_sum_.resize(n_forms);

		for(std::size_t i = 0; i < n_forms; ++i) {
			mat_buffer_[i] = std::make_shared< moonolith::SparseMatrix<double> >(comm_);
			local_element_matrices_sum_[i] = 0.;

			switch(assembler_->type(i)) {

				case LocalAssembler::MASTER_X_SLAVE:
				{
					mat_buffer_[i]->set_size(to_n_dofs_, from_n_dofs_);
					break;
				}

				case LocalAssembler::SLAVE_X_SLAVE:
				{
					mat_buffer_[i]->set_size(to_n_dofs_, to_n_dofs_);
					break;
				}

				case LocalAssembler::MASTER_X_MASTER:
				{
					mat_buffer_[i]->set_size(from_n_dofs_, from_n_dofs_);
					break;
				}

				case LocalAssembler::SLAVE_X_MASTER:
				{
					mat_buffer_[i]->set_size(from_n_dofs_, to_n_dofs_);
					break;
				}

				default:
				{
					assert(false);
					break;
				}
			}
		}
	}

	void PetrovGalerkinAssembler::finalize_form(
		std::size_t buffer_num,
		PetrovGalerkinAssembler::SparseMatrix &mat
		)
	{

		assert(buffer_num < mat_buffer_.size());

		libMesh::dof_id_type n_dofs_on_proc_trial = 0;
		libMesh::dof_id_type n_dofs_on_proc_test  = 0;

		switch(assembler_->type(buffer_num)) {

			case LocalAssembler::MASTER_X_SLAVE:
			{
				n_dofs_on_proc_trial = from_n_local_dofs_;
				n_dofs_on_proc_test  = to_n_local_dofs_;
				break;
			}

			case LocalAssembler::SLAVE_X_SLAVE:
			{
				n_dofs_on_proc_trial = to_n_local_dofs_;
				n_dofs_on_proc_test  = to_n_local_dofs_;
				break;
			}

			case LocalAssembler::MASTER_X_MASTER:
			{
				n_dofs_on_proc_trial = from_n_local_dofs_;
				n_dofs_on_proc_test  = from_n_local_dofs_;
				break;
			}

			case LocalAssembler::SLAVE_X_MASTER:
			{
				n_dofs_on_proc_trial = to_n_local_dofs_;
				n_dofs_on_proc_test  = from_n_local_dofs_;
				break;
			}

			default:
			{
				assert(false);
				break;
			}
		}

		local2global_->redistribute(comm_, n_dofs_on_proc_trial, n_dofs_on_proc_test, *mat_buffer_[buffer_num]);

		SizeType m_max_row_entries = mat_buffer_[buffer_num]->local_max_entries_x_col();
		comm_.all_reduce(&m_max_row_entries, 1, moonolith::MPIMax());

		SparseMatrix mat_x = utopia::local_sparse(n_dofs_on_proc_test, n_dofs_on_proc_trial, m_max_row_entries);

		{
			utopia::Write<SparseMatrix> write(mat_x);
			for (auto it = mat_buffer_[buffer_num]->iter(); it; ++it) {
				mat_x.set(it.row(), it.col(), *it);

			}
		}

		if(opts_.n_var == 1) {
			mat = std::move(mat_x);
			return;
		}

		auto s_mat_x = local_size(mat_x);
		mat = local_sparse(s_mat_x.get(0), s_mat_x.get(1), opts_.n_var * m_max_row_entries);

		utopia::Write<SparseMatrix> w_mat(mat);
		utopia::each_read(mat_x, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
			for(utopia::SizeType d = 0; d < opts_.n_var; ++d) {
				mat.set(i + d, j + d, value);
			}
		});
	}

	void PetrovGalerkinAssembler::print_stats()
	{
        auto qb_assembler = std::dynamic_pointer_cast<QuadratureBasedAssembler>(assembler_);
        if(qb_assembler) {
            double total_intersection_volume = qb_assembler->get_q_builder().get_total_intersection_volume();
            double volumes[2] = { local_element_matrices_sum_[0], total_intersection_volume };
            long isect_stats[2] = { n_intersections_, n_false_positives_};
            comm_.all_reduce(volumes, 2, moonolith::MPISum());
            comm_.all_reduce(isect_stats, 2, moonolith::MPISum());

            if(comm_.is_root()) {
                std::cout << "sum(B): "
                          << volumes[0]
                          << ", vol(I): "
                          << volumes[1]
                          << "\nn_candidates: "
                          << (isect_stats[0] + isect_stats[1])
                          << ", n_intersections: "
                          << isect_stats[0]
                          << ", n_false_positives: "
                          << isect_stats[1]
                          << std::endl;
            }

            // moonolith::root_describe("vol(I_local) : " + std::to_string(total_intersection_volume), comm, std::cout);
        }

	}

}
}
#ifndef UTOPIA_LIBMESH_ASSEMBLER_HPP
#define UTOPIA_LIBMESH_ASSEMBLER_HPP 

#include "utopia_libmesh_AssemblyContext.hpp"

namespace utopia {
	class LibMeshAssembler {
	public:
		typedef utopia::USparseMatrix GlobalMatrix;
		typedef utopia::UVector GlobalVector;
		typedef UTOPIA_SCALAR(GlobalVector) Scalar;

		LibMeshAssembler()
		: verbose_(Utopia::instance().verbose())
		{}

		template<class Expr>
		bool assemble(const Expr &expr, Scalar &val)
		{
			//perf
			Chrono c;
			c.start();

			typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
			typedef typename TraitsT::Matrix ElementMatrix;
			typedef typename TraitsT::Vector ElementVector;

			static const int Backend = TraitsT::Backend;

			const auto &space = find_space<LibMeshFunctionSpace>(expr);
			const auto &dof_map = space.dof_map();
			auto &m = space.mesh();


			val = 0.;

			for(auto it = elements_begin(m); it != elements_end(m); ++it) {
				init_context_on(expr, (*elements_begin(m))->id());

				if(it != elements_begin(m)) {
					reinit_context_on(expr, (*it)->id());
				}

				Number<Scalar> el_val = 0.;

				FormEvaluator<LIBMESH_TAG> eval;
				eval.eval(expr, el_val, ctx_);

				if(ctx_.has_assembled()) {
					val += el_val;
				}
			}

			m.comm().sum(val);

			//perf
			c.stop();

			if(verbose_) {
				std::cout << "assemble: value" << std::endl;
				std::cout << c << std::endl;
			}

			return false;
		}


		template<class Expr>
		bool assemble(const Expr &expr, GlobalMatrix &mat, GlobalVector &vec, const bool apply_constraints = false)
		{
			//perf
			Chrono c;
			c.start();

			typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
			typedef typename TraitsT::Matrix ElementMatrix;
			typedef typename TraitsT::Vector ElementVector;

			static const int Backend = TraitsT::Backend;

			const auto &space = find_space<LibMeshFunctionSpace>(expr);
			const auto &dof_map = space.dof_map();
			auto &m = space.mesh();

			auto s_m = size(mat);

			if(empty(mat) || s_m.get(0) != dof_map.n_dofs() || s_m.get(1) != dof_map.n_dofs()) {
				auto nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
					*std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));

				mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
			} else {
				mat *= 0.;
			}

			if(empty(vec) || size(vec).get(0) != dof_map.n_dofs()) {
				vec = local_zeros(dof_map.n_local_dofs());
			} else {
				vec.set(0.);
			}

			{
				Write<GlobalMatrix> w_m(mat);
				Write<GlobalVector> w_v(vec);

				ElementMatrix el_mat;
				ElementVector el_vec;

				init_context_on(expr, (*elements_begin(m))->id());

				for(auto it = elements_begin(m); it != elements_end(m); ++it) {
					if(it != elements_begin(m)) {
						reinit_context_on(expr, (*it)->id());
					}

					el_mat.implementation().zero();
					el_vec.implementation().zero();

					FormEvaluator<LIBMESH_TAG> eval;
					eval.eval(expr, el_mat, el_vec, ctx_);

					std::vector<libMesh::dof_id_type> dof_indices;
					dof_map.dof_indices(*it, dof_indices);

					if(ctx_.has_assembled()) {
						if(apply_constraints) {
							dof_map.heterogenously_constrain_element_matrix_and_vector(el_mat.implementation(), el_vec.implementation(), dof_indices);
						}

						add_matrix(el_mat.implementation(), dof_indices, dof_indices, mat);
						add_vector(el_vec.implementation(), dof_indices, vec);
					}
				}
			}

			//perf
			c.stop();
			if(verbose_) {
				std::cout << "assemble: lhs == rhs" << std::endl;
				std::cout << c << std::endl;
			}

			return true;
		}


		template<class Expr>
		bool assemble(const Expr &expr, GlobalMatrix &mat)
		{
			//perf
			Chrono c;
			c.start();

			typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
			typedef typename TraitsT::Matrix ElementMatrix;

			static const int Backend = TraitsT::Backend;

			const auto &space = find_space<LibMeshFunctionSpace>(expr);
			const auto &dof_map = space.dof_map();
			auto &m = space.mesh();

			auto s_m = size(mat);


			if(empty(mat) || s_m.get(0) != dof_map.n_dofs() || s_m.get(1) != dof_map.n_dofs()) {
				SizeType nnz_x_row = 0;
				if(!dof_map.get_n_nz().empty()) {
					// nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
					// 	*std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));

					nnz_x_row = 
						*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()) + 
						*std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end());
				}

				mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
			} else {
				mat *= 0.;
			}

			{
				Write<GlobalMatrix> w_m(mat, utopia::GLOBAL_ADD);

				if(elements_begin(m) != elements_end(m)) {

					ElementMatrix el_mat;
					init_context_on(expr, (*elements_begin(m))->id());

					for(auto it = elements_begin(m); it != elements_end(m); ++it) {
						if(it != elements_begin(m)) {
							reinit_context_on(expr, (*it)->id());
						}

						el_mat.implementation().zero();

						FormEvaluator<LIBMESH_TAG> eval;
						eval.eval(expr, el_mat, ctx_, true);

						std::vector<libMesh::dof_id_type> dof_indices;
						dof_map.dof_indices(*it, dof_indices);

						if(ctx_.has_assembled()) {
							add_matrix(el_mat.implementation(), dof_indices, dof_indices, mat);
						}
					}
				}
			}

			//perf
			c.stop();

			if(verbose_) {
				std::cout << "assemble: lhs" << std::endl;
				std::cout << c << std::endl;
			}

			return true;
		}


		template<class Expr>
		bool assemble(const Expr &expr, GlobalVector &vec, const bool apply_constraints = false)
		{

			//perf
			Chrono c;
			c.start();

			typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
			typedef typename TraitsT::Vector ElementVector;

			static const int Backend = TraitsT::Backend;

			const auto &space = find_space<LibMeshFunctionSpace>(expr);
			const auto &dof_map = space.dof_map();
			auto &m = space.mesh();

			if(empty(vec) || size(vec).get(0) != dof_map.n_dofs()) {
				// vec = local_zeros(dof_map.n_local_dofs());
				vec = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list()); 
			} else {
				vec *= 0.;
			}

			{
				Write<GlobalVector> w_v(vec, utopia::GLOBAL_ADD);
				ElementVector el_vec;

				if(elements_begin(m) != elements_end(m)) {
					init_context_on(expr, (*elements_begin(m))->id());

					for(auto it = elements_begin(m); it != elements_end(m); ++it) {
						if(it != elements_begin(m)) {
							reinit_context_on(expr, (*it)->id());
						}

						el_vec.implementation().zero();

						FormEvaluator<LIBMESH_TAG> eval;
						eval.eval(expr, el_vec, ctx_, true);

						std::vector<libMesh::dof_id_type> dof_indices;
						dof_map.dof_indices(*it, dof_indices);

						if(ctx_.has_assembled()) {
							add_vector(el_vec.implementation(), dof_indices, vec);
						}
					}
				}
			}


			//perf
			c.stop();

			if(verbose_) {
				std::cout << "assemble: rhs" << std::endl;
				std::cout << c << std::endl;
			}

			return true;
		}

		template<class Expr>
		void init_context_on(const Expr &expr, const int element_id)
		{	
			ctx_.set_current_element(element_id);
			ctx_.set_has_assembled(false);
			ctx_.init(expr);
		}

		template<class Expr>
		void reinit_context_on(const Expr &expr, const int element_id)
		{
			ctx_.set_current_element(element_id);
			ctx_.set_has_assembled(false);
			ctx_.reinit(expr);
		}

	private:
		AssemblyContext<LIBMESH_TAG> ctx_;
		bool verbose_;
	};
}

#endif //UTOPIA_LIBMESH_ASSEMBLER_HPP

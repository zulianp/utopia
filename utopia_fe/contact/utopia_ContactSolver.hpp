#ifndef UTOPIA_STEADY_CONTACTHPP
#define UTOPIA_STEADY_CONTACTHPP

#include "utopia.hpp"
#include "utopia_materials.hpp"
#include "utopia_Contact.hpp"
#include "utopia_Mechanics.hpp"
#include "utopia_SemiGeometricMultigrid.hpp"
#include "utopia_petsc_TaoSolver.hpp"
#include "utopia_ProjectedGradient.hpp"

#include "utopia_libmesh.hpp"

#include "libmesh/nemesis_io.h"
#include "libmesh/exodusII_io.h"

#include <memory>
#include <fstream>

namespace utopia {

	template<class Matrix, class Vector>
	class ContactSolver {
	public:
		DEF_UTOPIA_SCALAR(Matrix)
		typedef utopia::ProductFunctionSpace<LibMeshFunctionSpace> FunctionSpaceT;
		typedef libMesh::Nemesis_IO Exporter;
		// typedef libMesh::ExodusII_IO Exporter;

		ContactSolver(
			const std::shared_ptr<FunctionSpaceT> &V,
			const std::shared_ptr<ElasticMaterial<Matrix, Vector>> &material,
			const ContactParams &params)
		: V_(V),
		  material_(material),
		  params_(params),
		  first_(true),
		  tol_(1e-10),
		  debug_output_(false),
		  force_direct_solver_(false),
		  bypass_contact_(false),
		  exit_on_contact_solve_failure_(true),
		  sol_to_gap_on_contact_bdr_(false),
		  max_outer_loops_(20),
		  use_ssn_(false),
		  use_pg_(false),
		  max_non_linear_iterations_(30)
		{
			io_ = std::make_shared<Exporter>(V_->subspace(0).mesh());

			output_path_ = utopia::Utopia::instance().get("output_path");

			if(!output_path_.empty()) {
				output_path_ += "/";
			}

			output_path_ += "contact_sol.e";
			// linear_solver_ = std::make_shared<Factorization<Matrix, Vector>>();
			auto  iterative_solver = std::make_shared<BiCGStab<Matrix, Vector>>();
			// auto iterative_solver = std::make_shared<GaussSeidel<Matrix, Vector>>();

			iterative_solver->atol(1e-18);
			iterative_solver->stol(1e-17);
			iterative_solver->rtol(1e-8);
			iterative_solver->max_it(4000);
			linear_solver_ = iterative_solver;

			n_exports = 0;
 
			tao_.set_type("tron"); //REMOVED_TRILINOS
			tao_.set_ksp_types("bcgs", "jacobi", " "); //REMOVED_TRILINOS
		}

		void set_tol(const Scalar tol)
		{
			tol_ = tol;
		}

		void set_material(const std::shared_ptr< ElasticMaterial<Matrix, Vector> > &material)
		{
			material_ = material;
		}

		virtual ~ContactSolver() {}

		void update_contact(const Vector &x)
		{
			auto &V_0 = V_->subspace(0);

			if(bypass_contact_) {
				if(!contact_.initialized) {
					contact_.init_no_contact(
						utopia::make_ref(V_0.mesh()),
				    	utopia::make_ref(V_0.dof_map()));
				}

				return;
			}

			deform_mesh(V_0.mesh(), V_0.dof_map(), x);

			contact_.init(
				utopia::make_ref(V_0.mesh()),
				utopia::make_ref(V_0.dof_map()),
				params_
			);

			deform_mesh(V_0.mesh(), V_0.dof_map(), -x);

			auto mg = std::dynamic_pointer_cast<SemiGeometricMultigrid>(linear_solver_);

			if(mg) {
				mg->update_contact(contact_);
			}
		}

		bool solve_steady()
		{
			initialize();
			if(!solve_contact()) {
				assert(false);
				return false;
			}

			convert(x_, *V_->subspace(0).equation_system().solution);
			io_->write_equation_systems(output_path_, V_->subspace(0).equation_systems());

			finalize();
			return true;
		}

		bool solve_dynamic(const int n_time_steps)
		{
			initialize();

			n_exports = 0;



			convert(x_, *V_->subspace(0).equation_system().solution);
			io_->write_timestep(output_path_, V_->subspace(0).equation_systems(), n_exports + 1, n_exports);

			++n_exports;
			for(int t = 0; t < n_time_steps; ++t) {
				std::cout << "-------------------------------------"<< std::endl;
				std::cout << "time_step: " << t << std::endl;

				first_ = true;
				if(!solve_contact() && exit_on_contact_solve_failure_) return false;
				next_step();

				convert(x_, *V_->subspace(0).equation_system().solution);
				io_->write_timestep(output_path_, V_->subspace(0).equation_systems(), n_exports + 1, n_exports);



				++n_exports;
				std::cout << "-------------------------------------"<< std::endl;
			}

			finalize();
			return true;
		}

		bool solve_contact()
		{

			Vector old_sol = x_;

			for(int i = 0; i < max_outer_loops_; ++i) {
				contact_is_outdated_ = true;
				solve_contact_in_current_configuration();

				const double diff = norm2(old_sol - x_);

				if(debug_output_) {
					convert(x_, *V_->subspace(0).equation_system().solution);
					io_->write_timestep(output_path_, V_->subspace(0).equation_systems(), n_exports + 1, n_exports);
				}

				++n_exports;

				std::cout << "outer_loop: " << i << " diff: " << diff << std::endl;
				if(diff < tol_) {
					std::cout << "terminated at iteration " << i << " with diff " << diff << " <  " << tol_ << std::endl;
					break;
				} else {
					if(i + 1 == max_outer_loops_) {
						std::cerr << "[Warning] contact solver failed to converge with " << max_outer_loops_ << " loops under tolerance " << tol_ << std::endl;
						// assert(false);
						return false;
					}
				}

				old_sol = x_;
			}

			return true;
		}

		bool solve_contact_in_current_configuration()
		{
			bool converged = false;
			int iteration = 0;

			while(!converged) {

				if(!step()) return false;
				// if(material_->is_linear()) { break; }

				const double norm_inc = norm2(inc_c_);
				converged = norm_inc < tol_;

				std::cout << "iteration: " << iteration << " norm_inc: " << norm_inc << std::endl;
				++iteration;

				if(max_non_linear_iterations_ <= iteration) {
					std::cerr << "[Error] solver did not converge" << std::endl;
					return false;
				}
			}

			return true;
		}

		// virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient)
		virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient)
		{
			return material_->assemble_hessian_and_gradient(x, hessian, gradient);
		}

		// bool write_text(const std::string &path, const Matrix &mat)
		// {
		// 	int size = utopia::comm_size(mat);
		// 	int rank = utopia::comm_rank(mat);

		// 	int nnz = 0;
		// 	for(SizeType r = 0; r < size; ++r) {
		// 		if(r == 0) {
		// 			nnz = 0;
		// 			each_read(mat, [&nnz](const SizeType, const SizeType, const Scalar) {
		// 				++nnz;
		// 			});

		// 			MPI_Allreduce( MPI_IN_PLACE, &nnz, 1, MPI_INT, MPI_SUM, comm );
		// 		}

		// 		if(r == rank) {
		// 			std::ofstream os;

		// 			if(r == 0) {
		// 				os.open(path);
		// 				Size s = size(mat);
		// 				os << s.get(0) << " " << nnz << "\n";
		// 			} else {
		// 				os.open(path, std::ofstream::out | std::ofstream::app);
		// 			}

		// 			if(!os.good()) {
		// 				std::cerr << "invalid path: " << path << std::endl;
		// 				continue;
		// 			}

		// 			each_read(mat, [&os](const SizeType i, const SizeType j, const Scalar value) {
		// 				os << i << " " << j << " " << value << "\n";
		// 			});

		// 			os.flush();
		// 			os.close();
		// 		}

		// 		MPI_Barrier(comm);
		// 	}

		// 	return true;
		// }


		void qp_solve(Matrix &lhs, Vector &rhs, const BoxConstraints<Vector> &box_c, Vector &inc_c)
		{
			if(linear_solver_ && !contact_.has_contact()) {
				linear_solver_->solve(lhs, rhs, inc_c_);
				return;
			}

			if(sol_to_gap_on_contact_bdr_) {
				inc_c = e_mul(contact_.is_contact_node, *box_c.upper_bound());
			}

			// auto mg = std::dynamic_pointer_cast<SemiGeometricMultigrid>(linear_solver_);
			// if(!force_direct_solver_ && mg) {

			if(use_ssn_) {
			// 	SemismoothNewton<Matrix, Vector> newton(linear_solver_);
			// 	newton.verbose(true);
			// 	newton.max_it(40);
			// 	newton.set_box_constraints(box_c);
			// 	newton.solve(lhs, rhs, inc_c);
			// } else {
				// SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL> ssn;
				SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL> ssn(linear_solver_);  //REMOVED_TRILINOS
				// SemismoothNewton<Matrix, Vector> ssn(linear_solver_); 
				// ssn.verbose(true);
				ssn.max_it(40); //REMOVED_TRILINOS
				ssn.atol(1e-14); //REMOVED_TRILINOS
				ssn.rtol(1e-8); //REMOVED_TRILINOS
				ssn.stol(1e-8); //REMOVED_TRILINOS


				ssn.set_box_constraints(box_c); //REMOVED_TRILINOS
				ssn.solve(lhs, rhs, inc_c); //REMOVED_TRILINOS
			} else if(use_pg_) {

				ProjectedGradient<Matrix, Vector> pg;
				// SemismoothNewton<Matrix, Vector> pg(linear_solver_);
				// pg.verbose(true);
				pg.max_it(size(inc_c).get(0) * 20);
				pg.atol(1e-14);
				pg.rtol(1e-8);
				pg.stol(1e-8);

				pg.set_box_constraints(box_c);
				pg.solve(lhs, rhs, inc_c);

			} else {
				Chrono c;
				c.start();

				// static int n_stuff = 0;

				// if(n_stuff++ >= 2) {
				// 	std::cout << "Writing..." << std::flush;
					// write("rhs_" + std::to_string(size(lhs).get(0)) + ".m", rhs);
					// write("g_" + std::to_string(size(lhs).get(0)) + ".m", *box_c.upper_bound());
					// write_text("lhs_" + std::to_string(size(lhs).get(0)) + ".txt", lhs);
					// write("lhs_" + std::to_string(size(lhs).get(0)) + ".m", lhs);
				// 	std::cout << "done." << std::endl;
				// 	exit(0);
				// }

				QuadraticFunction<Matrix, Vector> fun(make_ref(lhs), make_ref(rhs));
				//(linear_solver_);
				tao_.set_box_constraints(box_c); //REMOVED_TRILINOS
				// tao_.set_ksp_types(KSPPREONLY, PCLU, "mumps");
				// tao_.set_type("gpcg");


				tao_.solve(fun, inc_c); //REMOVED_TRILINOS

				force_direct_solver_ = false;

				c.stop();

				std::cout << "Solve " << c << std::endl;
			}
			// }
		}

		bool step()
		{
			assert(x_.implementation().has_ghosts());
			x_.implementation().update_ghosts();

			if(contact_is_outdated_) {
				update_contact(x_);
				xc_ *= 0.;
				lagrange_multiplier_ *= 0.;
				contact_is_outdated_ = false;
			}

			if(!assemble_hessian_and_gradient(x_, H_, g_)) {
				assert(false);
				return false;
			}

			double norm_g = norm2(g_);
			std::cout << "norm_g: " << norm_g << std::endl;

			//handle transformations
			const auto &T = contact_.complete_transformation;

			gc_ = transpose(T) * g_;
			//change sign to negative gradient
			gc_ *= -1.;
			Hc_ = transpose(T) * H_ * T;

			apply_boundary_conditions(V_->subspace(0).dof_map(), Hc_, gc_);

			if(!first_) {
				apply_zero_boundary_conditions(V_->subspace(0).dof_map(), gc_);
			}

			inc_c_ *= 0.;
			qp_solve(Hc_, gc_, make_upper_bound_constraints(std::make_shared<Vector>(contact_.gap - xc_)), inc_c_);

			xc_ += inc_c_;
			x_ += T * inc_c_;

			first_ = false;
			return true;
		}

		void reset()
		{
			first_ = true;
			n_exports = 0;
		}

		virtual void initialize()
		{
			reset();
			auto &dof_map = V_->subspace(0).dof_map();
			x_ = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());
			inc_c_ = local_zeros(local_size(x_));
			xc_  = local_zeros(local_size(x_));
			lagrange_multiplier_ = local_zeros(local_size(x_));
		}

		virtual void finalize()
		{

		}

		const Contact &contact() const
		{
			return contact_;
		}

		ElasticMaterial<Matrix, Vector> &material()
		{
			return *material_;
		}


		FunctionSpaceT &space()
		{
			return *V_;
		}

		const FunctionSpaceT &space() const
		{
			return *V_;
		}

		const Vector &displacement() const
		{
			return x_;
		}

		Vector &displacement()
		{
			return x_;
		}

		void debug_output(const bool val)
		{
			debug_output_ = val;
		}


		virtual void next_step()
		{

		}


		inline const std::shared_ptr<ExternalForce> &external_force_fun() const
		{
			return external_force_fun_;
		}

		inline void set_external_force_fun(const std::shared_ptr<ExternalForce> &external_force_fun)
		{
			external_force_fun_ = external_force_fun;
		}

		void set_linear_solver(const std::shared_ptr<LinearSolver<Matrix, Vector> > &linear_solver)
		{
			linear_solver_ = linear_solver;
		}

		void set_bypass_contact(const bool val)
		{
			bypass_contact_ = val;
		}
		void set_max_outer_loops(const int val)
		{
			max_outer_loops_ = val;
		}

		void set_max_non_linear_iterations(const int val)
		{
			max_non_linear_iterations_ = val;
		}

		void set_use_ssn(const bool val)
		{
			use_ssn_ = val;
		}

		void set_use_pg(const bool val) {
			use_pg_ = val;
		}

		void set_exit_on_contact_solve_failure(const bool val)
		{
			exit_on_contact_solve_failure_ = val;
		}

		void set_sol_to_gap_on_contact_bdr(const bool val) {
			sol_to_gap_on_contact_bdr_ = val;
		}

		TaoSolver<Matrix, Vector> &tao() //REMOVED_TRILINOS
		{								 //REMOVED_TRILINOS
			return tao_;				 //REMOVED_TRILINOS
		}								//REMOVED_TRILINOS

		virtual bool stress(const Vector &x, Vector &result) {
			return material_->stress(x, result);
		}

	private:
		std::shared_ptr<FunctionSpaceT> V_;
		std::shared_ptr<ElasticMaterial<Matrix, Vector>> material_;
		std::shared_ptr<ExternalForce> external_force_fun_;
		ContactParams params_;
		bool first_;
		bool contact_is_outdated_;

		Scalar tol_;

		std::shared_ptr<LinearSolver<Matrix, Vector> > linear_solver_;


		Matrix H_;
		Vector g_;
		Vector x_;

		Matrix Hc_;
		Vector gc_;
		Vector inc_c_;
		Vector xc_;
		Vector rhs_;

		Contact contact_;

		Vector inactive_set_;
		Vector active_set_;
		USMatrix A_, I_;

		Vector lagrange_multiplier_;

		//additional vectors
		// USMatrix internal_mass_matrix_;

		std::shared_ptr<Exporter> io_;
		int n_exports;

		std::string output_path_;
		bool debug_output_;
		bool force_direct_solver_;
		bool bypass_contact_;
		bool exit_on_contact_solve_failure_;
		bool sol_to_gap_on_contact_bdr_;

		int max_outer_loops_;

		TaoSolver<Matrix, Vector> tao_; //REMOVED_TRILINOS
		bool use_ssn_, use_pg_;

		int max_non_linear_iterations_;
	};

	void run_steady_contact(libMesh::LibMeshInit &init);

}


#endif //UTOPIA_STEADY_CONTACTHPP
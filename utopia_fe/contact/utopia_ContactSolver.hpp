#ifndef UTOPIA_STEADY_CONTACTHPP
#define UTOPIA_STEADY_CONTACTHPP


#include "utopia_libmesh_Types.hpp"
#include "utopia_fe_config.hpp"

#ifndef WITH_TRILINOS_ALGEBRA

#include "utopia.hpp"
#include "utopia_materials.hpp"
#include "utopia_Contact.hpp"
#include "utopia_ContactAssembler.hpp"
#include "utopia_Mechanics.hpp"
#include "utopia_SemiGeometricMultigrid.hpp"
#include "utopia_petsc_TaoSolver.hpp"
#include "utopia_ProjectedGradient.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_ContactStress.hpp"
#include "utopia_polymorphic_QPSolver.hpp"
#include "utopia_ConvertContactAssembler.hpp"

#include "utopia_libmesh.hpp"

#include "libmesh/nemesis_io.h"
#include "libmesh/exodusII_io.h"

#include <memory>
#include <fstream>
#include <cassert>

namespace utopia {

    template<class Matrix, class Vector>
    class ContactSolver : public Configurable {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef utopia::ProductFunctionSpace<LibMeshFunctionSpace> FunctionSpaceT;
        typedef libMesh::Nemesis_IO Exporter;
        // using ContactT = utopia::ContactAssembler;
        using ContactT = utopia::IContact;
        using ContactStressT = utopia::ContactStress<ProductFunctionSpace<LibMeshFunctionSpace>, Matrix, Vector>;

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
          plot_gap_(false),
          max_outer_loops_(20),
          use_ssn_(false),
          use_pg_(false),
          max_non_linear_iterations_(30),
          export_results_(false),
          aux_system_num_(-1)
        {
            contact_ = utopia::make_unique<Contact>();

            io_ = std::make_shared<Exporter>(V_->subspace(0).mesh());

            output_path_ = utopia::Utopia::instance().get("output_path");

            if(!output_path_.empty()) {
                output_path_ += "/";
            }

            output_path_ += "contact_sol.e";
            auto  iterative_solver = std::make_shared<GMRES<Matrix, Vector>>("bjacobi");

            iterative_solver->atol(1e-18);
            iterative_solver->stol(1e-17);
            iterative_solver->rtol(1e-6);
            iterative_solver->max_it(4000);
            linear_solver_ = iterative_solver;

            n_exports = 0;

            auto tao = std::make_shared<TaoQPSolver<Matrix, Vector>>();
            tao->tao_type("tron");
            tao->set_linear_solver(std::make_shared<GMRES<Matrix, Vector>>("bjacobi"));
            qp_solver_ = tao;            
        }

        virtual void read(Input &is) override {

            std::string discretization = "legacy";
            is.get("discretization", discretization);

            if(discretization != "legacy") {
                // contact_ = utopia::make_unique<ContactAssembler>();
                contact_ = utopia::make_unique<ConvertContactAssembler>();
            }

            is.get("qp-solver", [this](Input &in) {
                this->qp_solver_ = std::make_shared<PolymorphicQPSolver<USparseMatrix, UVector>>();
                this->qp_solver_->read(in);
            });

            is.get("outer-loop-tol", tol_);
            is.get("max-outer-loops", max_outer_loops_);
            is.get("plot-gap", plot_gap_);
            is.get("contact", *contact_);
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
                if(!contact_->initialized()) {
                    contact_->init_no_contact(
                        utopia::make_ref(V_0.mesh()),
                        utopia::make_ref(V_0.dof_map()));
                }

                return;
            }

            deform_mesh(V_0.mesh(), V_0.dof_map(), x);

            contact_->assemble(
                utopia::make_ref(V_0.mesh()),
                utopia::make_ref(V_0.dof_map()),
                params_
            );

            if(plot_gap_) {
                static int n_out = 0;
                UVector gap = e_mul(contact_->is_contact_node(), contact_->gap());
                UVector is_contact = contact_->is_contact_node();
                UVector normals = contact_->normals();
                UVector is_glue_node = contact_->is_glue_node();

                write("gap" + std::to_string(n_out) + ".e",          V_->subspace(0), gap);
                write("is_contact" + std::to_string(n_out) + ".e",   V_->subspace(0), is_contact);
                write("normals" + std::to_string(n_out) + ".e",      V_->subspace(0), normals);
                

                if(!empty(is_glue_node)) {
                    write("is_glue_node" + std::to_string(n_out) + ".e", V_->subspace(0), is_glue_node);
                }

                
                // plot_gap_ = false;
                n_out++;
            }

            deform_mesh(V_0.mesh(), V_0.dof_map(), -x);

            auto mg = std::dynamic_pointer_cast<SemiGeometricMultigrid>(linear_solver_);
        }

        bool solve_steady()
        {
            initialize();
            first_ = true;
            if(!solve_contact()) {
                assert(false);
                return false;
            }

            if(export_results_) {

                // if(plot_gap_) {
                //     convert(contact_->gap(), *V_->subspace(0).equation_system().solution);
                // } else {
                convert(x_, *V_->subspace(0).equation_system().solution);
                // }

                create_aux_system();
                update_aux_system(x_);

                io_->write_equation_systems(output_path_, V_->subspace(0).equation_systems());
            }

            finalize();
            return true;
        }

        bool solve_dynamic(const int n_time_steps)
        {
            initialize();

            n_exports = 0;

            if(export_results_) {
                convert(x_, *V_->subspace(0).equation_system().solution);
                io_->write_timestep(output_path_, V_->subspace(0).equation_systems(), n_exports + 1, n_exports);
            }

            ++n_exports;
            for(int t = 0; t < n_time_steps; ++t) {
                std::cout << "-------------------------------------"<< std::endl;
                std::cout << "time_step: " << t << std::endl;

                first_ = true;
                if(!solve_contact() && exit_on_contact_solve_failure_) return false;
                next_step();

                if(export_results_) {
                    convert(x_, *V_->subspace(0).equation_system().solution);
                    io_->write_timestep(output_path_, V_->subspace(0).equation_systems(), n_exports + 1, n_exports);
                }

                ++n_exports;
                std::cout << "-------------------------------------"<< std::endl;
            }

            finalize();
            return true;
        }

        bool solve_contact()
        {

#ifdef WITH_PETSC
            UTOPIA_PETSC_MEMUSAGE();
#endif //WITH_PETSC

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

                    std::cout << "outer_loop: " << (i + 1) << "/" << max_outer_loops_ << " diff: " << diff << std::endl;
                    if(diff < tol_) {
                        std::cout << "terminated at iteration " << i << " with diff " << diff << " <  " << tol_ << std::endl;
                        break;
                    } else {
                        if(i + 1 == max_outer_loops_) {
                            std::cerr << "[Warning] contact solver failed to converge with " << max_outer_loops_ << " loops under tolerance " << tol_ << std::endl;
                            return true;
                        }
                    }

                    old_sol = x_;
                }
            }

#ifdef WITH_PETSC
            UTOPIA_PETSC_MEMUSAGE();
#endif //WITH_PETSC
            return true;
        }

        bool solve_contact_in_current_configuration()
        {
            bool converged = false;
            int iteration = 0;

            while(!converged) {

                if(!step()) return false;

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

        void qp_solve(Matrix &lhs, Vector &rhs, const BoxConstraints<Vector> &box_c, Vector &inc_c)
        {
            if(linear_solver_ && !contact_->has_contact()) {
                linear_solver_->solve(lhs, rhs, inc_c_);
                return;
            }

            // if(plot_gap_) {
            //     inc_c = e_mul(contact_->is_contact_node(), *box_c.upper_bound());
            // }

            Chrono c;
            c.start();

            qp_solver_->set_box_constraints(box_c);
            qp_solver_->solve(lhs, rhs, inc_c);

            c.stop();

            std::cout << "Solve " << c << std::endl;
        }

        bool step()
        {
            assert(x_.implementation().has_ghosts());
            synchronize(x_);

            if(contact_is_outdated_) {
#ifdef WITH_PETSC
                std::cout << "pre contact ";
                UTOPIA_PETSC_MEMUSAGE();
#endif //WITH_PETSC
                update_contact(x_);
                xc_ *= 0.;
                lagrange_multiplier_ *= 0.;
                contact_is_outdated_ = false;
#ifdef WITH_PETSC
                std::cout << "post contact ";
                UTOPIA_PETSC_MEMUSAGE();
#endif //WITH_PETSC
            }

            //maybe apply boundary conditions here???
            if(material_->is_linear()) {
                apply_boundary_conditions(V_->subspace(0).dof_map(), x_);
            }

            if(!assemble_hessian_and_gradient(x_, H_, g_)) {
                assert(false);
                return false;
            }

            double norm_g = norm2(g_);
            std::cout << "norm_g: " << norm_g << std::endl;

            //handle transformations
            contact_->couple(g_, gc_);

            //change sign to negative gradient
            gc_ *= -1.;
            contact_->couple(H_, Hc_);


            std::cout << "applying bc.... " << std::flush;
            apply_boundary_conditions(V_->subspace(0).dof_map(), Hc_, gc_);


            if(!first_ || material_->is_linear()) {
                apply_zero_boundary_conditions(V_->subspace(0).dof_map(), gc_);
            } 

            std::cout << "done" << std::endl;

            inc_c_ *= 0.;

            if(contact_->has_glue()) {  
                auto constr = make_upper_bound_constraints(std::make_shared<Vector>(contact_->gap() - xc_));
                constr.fill_empty_bounds();

                auto &u = *constr.upper_bound();
                auto &l = *constr.lower_bound();

                Write<UVector> wu(u), wl(l);
                each_read(contact_->is_glue_node(), [&constr, &l, &u](const SizeType i, const double val) {
                    if(val > 1e-8) {
                        u.set(i, 0.0);
                        l.set(i, 0.0);
                    }   
                });

                qp_solve(Hc_, gc_, constr, inc_c_);
            } else {
                qp_solve(Hc_, gc_, make_upper_bound_constraints(std::make_shared<Vector>(contact_->gap() - xc_)), inc_c_);
            }

            xc_ += inc_c_;

            UVector inc;
            contact_->uncouple(inc_c_, inc);
            x_ += inc;

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

        const ContactT &contact() const
        {
            assert(contact_);
            return *contact_;
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

        void export_results(const bool val)
        {
            export_results_ = val;
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

        void plot_gap(const bool val) {
            plot_gap_ = val;
        }

        virtual bool stress(const Vector &x, Vector &result) {
            return material_->stress(x, result);
        }


        void create_aux_system()
        {
            if(aux_system_num_ > 0) {
                std::cout << "aux system already exists" << std::endl;
                return;
            }

            auto &V0 = V_->subspace(0);
            auto &es = V0.equation_systems();

            const int dim  = es.get_mesh().mesh_dimension();

            auto &aux = es.add_system<libMesh::LinearImplicitSystem>("contact_aux");
            aux_system_num_ = aux.number();

            auto &dof_map_main = V0.dof_map();
            auto order = dof_map_main.variable_order(0);

            FunctionSpaceT W;
            W *= LibMeshFunctionSpace(aux, aux.add_variable("stress_x", order, libMesh::LAGRANGE));
            W *= LibMeshFunctionSpace(aux, aux.add_variable("stress_y", order, libMesh::LAGRANGE));

            if(dim > 2) {
                W *= LibMeshFunctionSpace(aux, aux.add_variable("stress_z", order, libMesh::LAGRANGE));
            }

            aux.init();

            // auto m_form = inner(trial(W), test(W)) * dX;
            // utopia::assemble(m_form, aux_mass_matrix_);


            // aux_inv_mass_matrix_ = utopia::make_unique<GMRES<USparseMatrix, UVector>>("bjacobi");
            // aux_inv_mass_matrix_->update(utopia::make_ref(aux_mass_matrix_));


            // aux_inv_mass_vector_ = 1./sum(aux_mass_matrix_, 1);

        }

        void update_aux_system(UVector &x)
        {
            // synchronize(x);
            
            UVector s, unscaled_s;
            if(contact_stress_) {
                contact_stress_->assemble(x, unscaled_s);
                // contact_->apply_orthogonal_trafo(s, unscaled_s);
            } else {
                stress(x, s);
                contact_->remove_mass(s, unscaled_s);
            }

            auto &V0 = V_->subspace(0);
            auto &es = V0.equation_systems();
            auto &aux = es.get_system<libMesh::LinearImplicitSystem>("contact_aux");

            utopia::convert(unscaled_s, *aux.solution);
            aux.solution->close();

            double max_s = utopia::max(unscaled_s);
            double min_s = utopia::min(unscaled_s);

            std::cout << "min_s: " << min_s << std::endl;
            std::cout << "max_s: " << max_s << std::endl;
        }

        inline void set_qp_solver(const std::shared_ptr<QPSolver<Matrix, Vector>> &qp_solver)
        {
            qp_solver_ = qp_solver;
        }

        inline QPSolver<Matrix, Vector> &qp_solver()
        {
            assert(qp_solver_);
            return *qp_solver_;
        }

        inline const QPSolver<Matrix, Vector> &qp_solver() const
        {
            assert(qp_solver_);
            return *qp_solver_;
        }

        inline void set_contact_stress(const std::shared_ptr<ContactStressT> &stress)
        {
            contact_stress_ = stress;
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

        std::unique_ptr<IContact> contact_;

        Vector inactive_set_;
        Vector active_set_;
        USparseMatrix A_, I_;

        Vector lagrange_multiplier_;

        std::shared_ptr<Exporter> io_;
        int n_exports;

        std::string output_path_;
        bool debug_output_;
        bool force_direct_solver_;
        bool bypass_contact_;
        bool exit_on_contact_solve_failure_;
        bool plot_gap_;

        int max_outer_loops_;

        bool use_ssn_, use_pg_;

        int max_non_linear_iterations_;
        bool export_results_;

        int aux_system_num_;
        USparseMatrix aux_mass_matrix_;
        UVector 	  aux_inv_mass_vector_;
        std::unique_ptr<LinearSolver<USparseMatrix, UVector>> aux_inv_mass_matrix_;

        std::shared_ptr<QPSolver<Matrix, Vector>> qp_solver_;
        std::shared_ptr<ContactStressT> contact_stress_;
    };

    void run_steady_contact(libMesh::LibMeshInit &init);

}

#endif //WITH_TRILINOS_ALGEBRA
#endif //UTOPIA_STEADY_CONTACTHPP
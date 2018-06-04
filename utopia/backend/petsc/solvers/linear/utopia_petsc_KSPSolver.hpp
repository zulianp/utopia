#ifndef UTOPIA_PETSC_KSP_HPP
#define UTOPIA_PETSC_KSP_HPP

#include "utopia_Preconditioner.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_petsc_KSPWrapper.hpp"

#include <algorithm>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsys.h>

namespace utopia {
    
    PetscErrorCode UtopiaPCApplyShell(PC pc, Vec x, Vec y);
    PetscErrorCode MyKSPMonitor(KSP,PetscInt,PetscReal,void*);
    
    class KSPLog {
    public:
        Vec x_k_1;
        Vec x_k_2;
        
        KSPLog()
        : x_k_1(nullptr), x_k_2(nullptr)
        {}
        
        inline bool initialized() const
        {
            return x_k_1 != nullptr;
        }
        
        void init_from(Vec x)
        {
            destroy();
            
            VecDuplicate(x, &x_k_1);
            VecDuplicate(x, &x_k_2);
        }
        
        ~KSPLog()
        {
            destroy();
        }
        
        inline void destroy()
        {
            if(x_k_1) {
                VecDestroy(&x_k_1);
                x_k_1 = nullptr;
            }
            
            if(x_k_2) {
                VecDestroy(&x_k_2);
                x_k_2 = nullptr;
            }
        }
    };
    
    /**@ingroup     Linear
     * @brief       Class provides interface to Petsc KSP solvers \n
     *              For setting up basic parameters, one can use classic Petsc runtime options, e.g.
     *              To see all possibilities, please refer to:
     *                                                   * <a href="http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html">preconditioner types</a>
     *                                                   * <a href="http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType">solver types</a>
     *
     *              Setting own/utopia preconditioner, can be done as following:
     *              \snippet tests/utopia_SolverTest.cpp PetscKSPSolver solve example1
     *              Detailed information about preconditioners, can be found in  \ref precondotioners.
     */
    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class KSPSolver {};
    
    template<typename Matrix, typename Vector>
    class KSPSolver<Matrix, Vector, PETSC> : public PreconditionedSolver<Matrix, Vector>, public Smoother<Matrix, Vector> {
    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Preconditioner<Vector> Preconditioner;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::PreconditionedSolver<Matrix, Vector> PreconditionedSolver;
        typedef utopia::Smoother<Matrix, Vector> Smoother;
        
        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "only works with petsc types");
        
        KSPSolver(const Parameters params = Parameters())
        : ksp_(std::make_shared<KSPWrapper<Matrix, Vector>>(PETSC_COMM_WORLD, params))
        {
            ksp_type("bcgs");
            pc_type("jacobi");
            ksp_->set_initial_guess_non_zero(true);
            set_parameters(params);
            
        }
        
        KSPSolver(const std::shared_ptr<KSPWrapper<Matrix, Vector>> &w)
        : ksp_(w)
        {}
        
        inline void wrap(KSP &ksp)
        {
            ksp_ = std::make_shared< KSPWrapper<Matrix, Vector> >(ksp, false);
        }
        
        virtual ~KSPSolver() {}
        
        /**
         * @brief      Sets the parameters.
         *
         * @param[in]  params  The parameters
         */
        virtual void set_parameters(const Parameters params) override
        {
            PreconditionedSolver::set_parameters(params);
            ksp_->set_parameters(params);
        }
        
        /* @brief      Sets the choice of direct solver.
         *             Please note, in petsc, direct solver is used as preconditioner alone, with proper settings.
         *
         * @param[in]  pc_type  The type of direct solver.
         */
        void pc_type(const std::string &pc_type)
        {
            ksp_->pc_type(pc_type);
        }
        
        /**
         * @brief      Sets KSP type
         */
        void ksp_type(const std::string & ksp_type)
        {
            ksp_->ksp_type(ksp_type);
        }
        
        /**
         * @brief      Sets solver package for choice of direct solver.
         *
         * @param[in]  package  The solver package.
         */
        void solver_package(const std::string &package)
        {
            ksp_->solver_package(package);
        }
        
        /**
         * @brief      Returns type of direct solver.
         */
        std::string pc_type() const { return this->ksp_->pc_type();}
        
        /**
         * @brief      Returns ksp package type
         */
        std::string ksp_type() const { return this->ksp_->ksp_type();}
        
        /**
         * @brief      Returns type of solver package.
         */
        std::string solver_package() const { return this->ksp_->solver_package();}
        
        /**
         * @brief        Solve function for all PETS KSP solvers.
         *              It is also compatible with our own Utopia preconditioners.
         */
        bool apply(const Vector &b, Vector &x) override
        {
            ksp_->set_tolerances(this->rtol(), this->atol(), PETSC_DEFAULT, this->max_it());

            
            // is this proper place to do so??? 
            // this->set_ksp_options(ksp_->implementation());
            return ksp_->apply(b, x);
        }
        
        bool smooth(const Vector &rhs, Vector &x) override
        {
            return ksp_->smooth(this->sweeps(), rhs, x);
        }
        
        inline bool must_compute_cond_number() const
        {
            static const std::vector<std::string> types = {"bcgs", "cg", "gmres"};
            return (this->verbose() && in_array(ksp_->ksp_type(), types));
        }
        
        /**
         * @brief      Sets the default options for PETSC KSP solver. \n
         *             Default: BiCGstab
         *
         * @param      ksp   The ksp
         */
        virtual void set_ksp_options(KSP &ksp)
        {
            ksp_->copy_settings_to(ksp);
            set_monitor_options(ksp);   
        }
        
        virtual void attach_preconditioner(KSP &ksp) const {
            KSPWrapper<Matrix, Vector> w(ksp, false);
            auto delegate_ptr = std::dynamic_pointer_cast<DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());
            
            if(delegate_ptr) {
              if(ksp_->has_shell_pc()) {
                  m_utopia_warning_once("set_preconditioner sets jacobi if a delegate precond has been set and type is matshell");
              }
            } else if(this->get_preconditioner()) {
                auto shell_ptr = this->get_preconditioner().get();
                w.attach_shell_preconditioner(UtopiaPCApplyShell,
                                              shell_ptr,
                                              nullptr,
                                              nullptr
                                              );
            }
        }
        
        virtual void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override
        {
            PreconditionedSolver::set_preconditioner(precond);
            
            auto delegate_ptr = std::dynamic_pointer_cast<DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());
            
            if(delegate_ptr) {
                if(ksp_->has_shell_pc()) {
                    m_utopia_warning_once("set_preconditioner sets jacobi if a delegate precond has been set and type is matshell");
                    ksp_->pc_type("jacobi");
                }

            } else if(this->get_preconditioner()) {
                auto shell_ptr = this->get_preconditioner().get();
                ksp_->attach_shell_preconditioner(UtopiaPCApplyShell,
                                                  shell_ptr,
                                                  nullptr,
                                                  nullptr
                                                  );
            }
        }
        
        virtual void set_monitor_options(KSP &ksp) const
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            if(this->verbose()) {
                ierr = KSPMonitorSet(ksp,
                                     MyKSPMonitor,
                                     new KSPLog(),
                                     [](void **arg) -> PetscErrorCode {
                                         auto ksp_log = (KSPLog **)(arg);
                                         delete *ksp_log;
                                         *ksp_log = nullptr;
                                         return 0;
                                     });  assert(ierr == 0);
            }
            
            if(must_compute_cond_number()) {
                ierr = KSPSetComputeSingularValues(ksp, PETSC_TRUE); assert(ierr == 0);
            }
        }
        
        void handle_reset(const Matrix &op)
        {
            bool must_reset = true;
            // bool must_reset = false;
            
            // if(this->has_operator()) {
            //     auto s = size(*this->get_operator());
            //     auto other_s = size(op);
            //     if(s != other_s) {
            //         must_reset = true;
            //     }
            // }
            
            // if(op.implementation().communicator() != ksp_->communicator())
            // {
            //     must_reset = true;
            // }
            
            if(must_reset) {
                auto temp_ksp = std::make_shared<KSPWrapper<Matrix, Vector>>(op.implementation().communicator());
                temp_ksp->copy_settings_from(*ksp_);
                ksp_ = temp_ksp;

                if(this->get_preconditioner()) {
                  set_preconditioner(this->get_preconditioner());
                }
            }
        }
        
        virtual void update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec) override
        {
            handle_reset(*op);
            set_monitor_options(ksp_->implementation());
            

            PreconditionedSolver::update(op, prec);
            ksp_->update(*op, *prec);
        }
        
        virtual void update(const std::shared_ptr<const Matrix> &op) override
        {
            handle_reset(*op);
            

            PreconditionedSolver::update(op);
            set_monitor_options(ksp_->implementation());
            
            bool skip_set_operators = false;
            if(this->get_preconditioner()) {
                auto mat_prec = std::dynamic_pointer_cast< DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());
                if(mat_prec) {
                    ksp_->update(*op, *mat_prec->get_matrix());
                    skip_set_operators = true;
                }
            }
            
            if(!skip_set_operators) {
                ksp_->update(*op);
            }
        }
        
        inline KSPSolver &operator=(const KSPSolver &other)
        {
            if(this == &other) return *this;
            PreconditionedSolver::operator=(other);
            Smoother::operator=(other);
            
            ksp_ = std::make_shared<KSPWrapper<Matrix, Vector>>(other.ksp_->communicator());
            ksp_->copy_settings_from(*other.ksp_);
            return *this;
        }
        
        inline KSPSolver &operator=(KSPSolver &&other)
        {
            if(this == &other) return *this;
            PreconditionedSolver::operator=(std::move(other));
            Smoother::operator=(std::move(other));
            ksp_ = std::move(other.ksp_);
            return *this;
        }
        
        KSPSolver(const KSPSolver &other):
        PreconditionedSolver(other),
        Smoother(other),
        ksp_(std::make_shared<KSPWrapper<Matrix, Vector>>(other.ksp_->communicator()))
        {
            ksp_->copy_settings_from(*other.ksp_);
        }
        
        KSPSolver(KSPSolver &&other)
        : PreconditionedSolver(std::move(other)),
        Smoother(std::move(other)),
        ksp_(std::move(other.ksp_))
        {}
        
        virtual KSPSolver * clone() const override
        {
            return new KSPSolver(*this);
        }
        
        KSPWrapper<Matrix, Vector> &ksp()
        {
            return *ksp_;
        }
        
        const KSPWrapper<Matrix, Vector> &ksp() const
        {
            return *ksp_;
        }
        
    protected:
        std::shared_ptr<KSPWrapper<Matrix, Vector>> ksp_;
    };
    
}

#endif //UTOPIA_PETSC_KSP_HPP

#ifndef UTOPIA_PETSC_MULTIGRID_HPP
#define UTOPIA_PETSC_MULTIGRID_HPP

#include "utopia_Multigrid.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_petsc.hpp"

namespace utopia {


    template<class Matrix, class Vector>
    class Multigrid<Matrix, Vector, PETSC_EXPERIMENTAL> : public IterativeSolver<Matrix, Vector>, public LinearMultiLevel<Matrix, Vector> {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::LinearSolver<Matrix, Vector>        Solver;
        typedef utopia::IterativeSolver<Matrix, Vector>     IterativeSolver;
        typedef utopia::Smoother<Matrix, Vector>            Smoother;
        typedef utopia::Level<Matrix, Vector>               Level;
        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef utopia::MatrixTransfer<Matrix, Vector>      MatrixTransfer;
    public:
        void update(const std::shared_ptr<const Matrix> &op) override
        {
            if(!ksp_) {
                init_ksp(op->comm().get());
            }

            this->galerkin_assembly(op);
            KSPSetOperators(*ksp_, raw_type(*op), raw_type(*op));

            PC pc;
            KSPGetPC(*ksp_, &pc);

            std::size_t n_smoothers = this->n_levels()-1;
            for(std::size_t i = 0; i < n_smoothers; i++)
            {
                KSP smoother;
                PCMGGetSmoother(pc, i, &smoother);
                if(block_size_ > 1) {
                    const_cast<Matrix &>(this->level(i).A()).convert_to_mat_baij(block_size_);
                }

                KSPSetOperators(smoother, raw_type(this->level(i).A()), raw_type(this->level(i).A()));
            }
        }

        virtual bool apply(const Vector &rhs, Vector &x) override
        {
            if(this->verbose())
                this->init_solver("utopia/petsc Multigrid",  {" it.", "|| Au - b||"});

            KSPSolve(*ksp_, raw_type(rhs), raw_type(x));

            KSPConvergedReason  reason;
            PetscInt            its;
            KSPGetConvergedReason(*ksp_, &reason);
            KSPGetIterationNumber(*ksp_, &its);

            if(this->verbose())
                this->exit_solver(its, reason);
            //FIXME
            return true;
        }

        Multigrid(const std::shared_ptr<Smoother> &smoother    = nullptr,
                  const std::shared_ptr<Solver> &linear_solver = nullptr)
        : smoother_(smoother), linear_solver_(linear_solver), default_ksp_type_(KSPRICHARDSON), default_pc_type_(PCSOR), block_size_(1)
        {

        }

        inline void set_default_ksp_type(const KSPType &ksp_type)
        {
            default_ksp_type_ = ksp_type;
        }

        inline void set_default_pc_type(const PCType &pc_type)
        {
            default_pc_type_ = pc_type;
        }


        virtual void update_transfer(const SizeType level, const std::shared_ptr<Transfer> &t) override
        {
            LinearMultiLevel<Matrix, Vector>::update_transfer(level, t);
            aux_update_transfer(level);
        }


        inline void block_size(const PetscInt size)
        {
            block_size_ = size;
        }

        Multigrid * clone() const override
        {
            return new Multigrid();
        }


        void read(Input &in) override
        {
          LinearMultiLevel<Matrix, Vector>::read(in);
          IterativeSolver::read(in);

          in.get("block_size", block_size_);

          if(smoother_) {
              in.get("smoother", *smoother_);
          }
          if(linear_solver_) {
              in.get("coarse_solver", *linear_solver_);
          }

        }

        void print_usage(std::ostream &os) const override
        {
          LinearMultiLevel<Matrix, Vector>::print_usage(os);
          IterativeSolver::print_usage(os);

          this->print_param_usage(os, "block_size", "int", "Block size for systems.", "1");
          this->print_param_usage(os, "smoother", "Smoother", "Input parameters for all smoothers.", "-");
          this->print_param_usage(os, "coarse_solver", "LinearSolver", "Input parameters for coarse solver.", "-");
        }

    private:
        std::shared_ptr<Smoother> smoother_;
        std::shared_ptr<Solver>   linear_solver_;

        std::shared_ptr<KSP> ksp_;
        std::vector<Level> levels_;

        KSPType default_ksp_type_;
        PCType default_pc_type_;

        PetscInt block_size_;

        void aux_update_transfer(const SizeType level)
        {
            auto mat_transfer = dynamic_cast<MatrixTransfer *>(&this->transfer(level));

            if(!mat_transfer) {
                std::cerr << "[Error] FIXME incompatible transfer type" << std::endl;
                assert(false);
                return;
            }

            if(!ksp_) {
                init_ksp(mat_transfer->I().comm().get());
            } else {
                PC pc;
                KSPGetPC(*ksp_, &pc);
                Mat I = raw_type(mat_transfer->I());
                PCMGSetInterpolation(pc, level + 1, I);

                Mat R = raw_type(mat_transfer->R());
                PCMGSetRestriction(pc, level + 1, R);

                PCSetUp(pc);
            }
        }

        template<class Solver>
        bool set_solver(Solver &solver, KSP &ksp)
        {
            auto * casted = dynamic_cast< KSPSolver<Matrix, Vector> *>(&solver);

            if(casted) {
                // if(casted->ksp_type() == KSPCG || casted->ksp_type() == KSPFCG) {
                // 	KSPSetType(ksp, KSPFCG);
                // } else
                if(casted->ksp_type() == KSPGMRES || casted->ksp_type() == KSPFGMRES) {
                    KSPSetType(ksp, KSPFGMRES);
                } else {
                    if(casted->ksp_type() != default_ksp_type_) {
                        std::cerr << "[Warning] " << casted->ksp_type() << " not supported by petsc PCMG using fallback option KSPRICHARDSON" << std::endl;
                    }

                    KSPSetType(ksp, default_ksp_type_);
                }

                PC pc;
                KSPGetPC(ksp, &pc);

                if(!casted->get_preconditioner()) {
                    PCSetType(pc, casted->pc_type().c_str());
                } else {
                // 	//FIXME use PCShell
                    PCSetType(pc, default_pc_type_);
                }

                {
                    PCType pc_type;
                    PCGetType(pc, &pc_type);

                    KSPType ksp_type;
                    KSPGetType(ksp, &ksp_type);

                    std::cout << ksp_type << ", " << pc_type << std::endl;
                }

                // return false;
                return true;

            } else {
                auto * factor = dynamic_cast< Factorization<Matrix, Vector, PETSC> *>(&solver);

                if(factor) {
                    factor->strategy().set_ksp_options(ksp);
                    KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
                    return true;
                }

                // else {
                // 	return false;
                // }
            }

            return false;
        }

        void init_ksp(MPI_Comm comm)
        {
            ksp_ = std::shared_ptr<KSP>(new KSP, [](KSP *&ksp) { KSPDestroy(ksp); delete ksp; ksp = nullptr; });

            KSPCreate(comm, ksp_.get());
            KSPSetFromOptions(*ksp_);

            PC pc;
            KSPGetPC(*ksp_, &pc);
            PCSetType(pc, PCMG);

            PCMGSetLevels(pc, this->n_levels(), nullptr);
            // PCMGSetGalerkin(pc, PETSC_TRUE);
#if UTOPIA_PETSC_VERSION_LESS_THAN(3,8,0)
            PCMGSetGalerkin(pc, PETSC_FALSE);
#else
            PCMGSetGalerkin(pc, PC_MG_GALERKIN_NONE);
#endif
            KSPSetInitialGuessNonzero(*ksp_, PETSC_TRUE);

            const std::size_t n_smoothers = this->n_levels()-1;
            for (std::size_t i = 0; i < n_smoothers; i++)
            {
                KSP smoother;
                PC sm;
                PCMGGetSmoother(pc, i + 1, &smoother);

                bool user_solver = false;

                if(smoother_) {
                    user_solver = set_solver(*smoother_, smoother);
                    if(!user_solver) {
                        std::cerr << "[Warning] not using user smoother" << std::endl;
                    }
                }

                if(!user_solver) {
                    // KSPFGMRES, KSPGCG, or KSPRICHARDSON
                    KSPSetType(smoother, default_ksp_type_);
                    KSPGetPC(smoother, &sm);
                    PCSetType(sm, default_pc_type_);

                    // KSPSetInitialGuessNonzero(smoother, PETSC_TRUE);
                }

                auto mat_transfer = dynamic_cast<MatrixTransfer *>(&this->transfer(i));

                if(!mat_transfer) {
                    std::cerr << "[Error] FIXME incompatible transfer type" << std::endl;
                    assert(false);
                    return;
                }

                Mat I = raw_type(mat_transfer->I());
                PCMGSetInterpolation(pc, i + 1, I);

                Mat R = raw_type(mat_transfer->R());
                PCMGSetRestriction(pc, i + 1, R);

                if(linear_solver_) {
                    KSP coarse_solver;
                    PCMGGetCoarseSolve(pc, &coarse_solver);
                    bool user_solver = set_solver(*linear_solver_, coarse_solver);

                    if(!user_solver) {
                        std::cerr << "[Warning] not using user linear_solver" << std::endl;
                    }
                }
            }

#if UTOPIA_PETSC_VERSION_LESS_THAN(3,9,0)
            PCMGSetNumberSmoothUp(pc,   this->post_smoothing_steps());
            PCMGSetNumberSmoothDown(pc, this->pre_smoothing_steps());
#else
            m_utopia_warning_once("PCMGSetNumberSmooth{Up,Down} not available in petsc 3.9.0 find equivalent");
#endif

            KSPSetTolerances(*ksp_, this->rtol(), this->atol(), PETSC_DEFAULT, this->max_it());

            if(this->verbose()) {
                KSPMonitorSet(
                *ksp_,
                [](KSP, PetscInt iter, PetscReal res, void*) -> PetscErrorCode {
                    PrintInfo::print_iter_status({static_cast<Scalar>(iter), res});
                    return 0;
                },
                nullptr,
                nullptr);
            }
        }
    };

}

#endif //UTOPIA_PETSC_MULTIGRID_HPP

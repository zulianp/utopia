#ifndef UTOPIA_TR_SUBPROBLEM_KSP_TR_HPP
#define UTOPIA_TR_SUBPROBLEM_KSP_TR_HPP
#include "utopia_TRSubproblem.hpp"
#include <string>
#include <cstring>

namespace utopia
{
    /**
     * @brief      Interface to use petsc KSP in TR context \n
     *				There are 5 different possibilities, grouped as follows: \n
     *				1. gltr, stcg, nash, qcg: \n
                                             Uses preconditioned conjugate gradient to compute an approximate minimizer of the quadratic function.  \n
                                            \f$ q(s) = g^T * s + .5 * s^T * H * s   \f$  \n
                                            subject to  \f$ || s || <= delta  \f$ \n
                    2. cgne				: \n
                                            Applies the preconditioned conjugate gradient method to the normal equations without explicitly forming A^t*A. \n
                                              Therefore, it's very suitable to use in combination with TR_NormalEquation, LeastSquaresFunction and solution method: dogleg \n
                                              Note, that CGNE is a general-purpose non-symmetric method. \n
     */
    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class KSP_TR {};


    template<typename Matrix, typename Vector>
    class KSP_TR<Matrix, Vector, PETSC>: public TRSubproblem<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef utopia::KSPSolver<Matrix, Vector>           KSPSolver;
        typedef utopia::TRSubproblem<Matrix, Vector>        TRSubproblem;

        typedef utopia::Preconditioner<Vector>              Preconditioner;
        typedef utopia::IterativeSolver<Matrix, Vector>     IterativeSolver;

        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "utopia::KSP_TR:: only works with petsc types");

    public:

        KSP_TR(const bool redundant_flg = false): TRSubproblem(), redundant_solve_flg_(redundant_flg)//,{"stcg", "nash", "cgne", "gltr", "qcg"})
        {
            this->ksp_type("stcg");
            this->pc_type("jacobi");
        }

        KSP_TR(const std::string type, const std::string pc_type = "jacobi", const bool redundant_flg = false): TRSubproblem(), redundant_solve_flg_(redundant_flg)
        {
            this->pc_type(pc_type);
            this->ksp_type(type);
        }


        virtual ~KSP_TR(){}

        virtual void read(Input &in) override
        {
            TRSubproblem::read(in);

            std::string pc_type_aux;
            std::string ksp_type_aux;

            in.get("pc_type", pc_type_aux);
            in.get("ksp_type", ksp_type_aux);
            in.get("redundant_solve_flg", redundant_solve_flg_);

            ksp_.pc_type(pc_type_aux);
            ksp_.ksp_type(ksp_type_aux);

        }

        virtual void ksp_type(const std::string & ksp_type_name)
        {
            if(!redundant_solve_flg_)
            {
                ksp_.ksp_type(ksp_type_name);
            }
            else
            {
                ksp_.ksp_type("preonly");
                ksp_.pc_type("redundant");

                // setting up inner solver
                PC pc_redundant;
                KSPGetPC(ksp_.implementation(), &pc_redundant);

                KSP innerksp;
                PCRedundantGetKSP(pc_redundant, &innerksp);
                KSPSetType(innerksp, ksp_type_name.c_str());

            }
        }

        virtual void pc_type(const std::string & pc_type_name)
        {
            if(!redundant_solve_flg_)
            {
                ksp_.pc_type(pc_type_name);
            }
            else
            {
                ksp_.ksp_type("preonly");
                ksp_.pc_type("redundant");

                // setting up inner solver
                PC pc_redundant;
                KSPGetPC(ksp_.implementation(), &pc_redundant);

                KSP innerksp; PC inner_pc;
                PCRedundantGetKSP(pc_redundant, &innerksp);
                KSPGetPC(innerksp, &inner_pc);
                PCSetType(inner_pc, pc_type_name.c_str());
            }
        }


        virtual void print_usage(std::ostream &os) const override
        {
            TRSubproblem::print_usage(os);

            this->print_param_usage(os, "pc_type", "string", "Type of petsc preconditioner.", "jacobi");
            this->print_param_usage(os, "ksp_type", "string", "Type of ksp solver to be used.", "stcg");
        }

        void number_of_parallel_solves(const SizeType & number)
        {
            if(redundant_solve_flg_)
            {
                PC pc_redundant;
                KSPGetPC(ksp_.implementation(), &pc_redundant);

                PetscBool flg;
                PetscObjectTypeCompare((PetscObject)pc_redundant,PCREDUNDANT,&flg);

                if(!flg)
                {
                    PCSetType(pc_redundant, PCREDUNDANT);
                    PCRedundantSetNumber(pc_redundant, number);

                    this->ksp_type("stcg");
                    this->pc_type("jacobi");
                }
                else
                {
                    PCRedundantSetNumber(pc_redundant, number);
                }
            }
            else
            {
                utopia_error("KSP_TR::number_of_parallel_solves:: If you want to specify number of solves, use redundant solver first.");
            }
        }

        bool redundant() const
        {
            return redundant_solve_flg_;
        }


    public:
        virtual KSP_TR<Matrix, Vector, PETSC>* clone() const override {
            return new KSP_TR<Matrix, Vector, PETSC>(*this);
        }

        virtual bool apply(const Vector &b, Vector &x) override
        {
            ksp_.apply(b, x);
            return true;
        }

        /**
         * @brief      Update function.
         */
        virtual void update(const std::shared_ptr<const Matrix> &op) override
        {
            ksp_.update(op);
            set_ksp_options();
        }


        virtual void set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
        {
            if(!redundant_solve_flg_)
            {
                ksp_.set_preconditioner(precond);
            }
            else
            {
                utopia_error("KSP_TR::set_preconditioner not supported for redundant solver.");
            }
        }

        // necessary, as ksp_ is resetting options with every apply...
        virtual void atol(const Scalar & atol_in) override
        {
            TRSubproblem::atol(atol_in);
            ksp_.atol(atol_in);
        }

        virtual void stol(const Scalar & stol_in)  override
        {
            TRSubproblem::stol(stol_in);
            ksp_.stol(stol_in);
        }

        virtual void rtol(const Scalar & rtol_in) override
        {
            TRSubproblem::rtol(rtol_in);
            ksp_.rtol(rtol_in);
        }

        virtual void max_it(const SizeType & max_it_in) override
        {
            TRSubproblem::max_it(max_it_in);
            ksp_.max_it(max_it_in);
        }

        virtual void verbose(const bool & verbose_in) override
        {
            TRSubproblem::verbose(verbose_in);
            ksp_.verbose(verbose_in);
        }


    private:

           /**
         * @brief      Sets the default options for PETSC KSP solver. \n
         *             Default: ST-CG
         *
         * @param      ksp   The ksp
         */
        void set_ksp_options()
        {
            if(!redundant_solve_flg_)
            {
                PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
                ierr = KSPSetFromOptions(ksp_.implementation());

                ierr = KSPSetType(ksp_.implementation(), ksp_.ksp_type().c_str());
                ierr = KSPSetInitialGuessNonzero(ksp_.implementation(), PETSC_TRUE);

                if(!ksp_.get_preconditioner())
                {
                    PC pc;
                    ierr = KSPGetPC(ksp_.implementation(), &pc);
                    ierr = PCSetType(pc, ksp_.pc_type().c_str());
                }

    #if UTOPIA_PETSC_VERSION_LESS_THAN(3,8,0)
                if(ksp_.ksp_type() == "qcg")
                    ierr = KSPQCGSetTrustRegionRadius(ksp_.implementation(), this->current_radius());
                else if(ksp_.ksp_type() == "gltr")
                    ierr = KSPGLTRSetRadius(ksp_.implementation(), this->current_radius());
                else if(ksp_.ksp_type() == "nash")
                    ierr = KSPNASHSetRadius(ksp_.implementation(), this->current_radius());
                else
                    ierr = KSPSTCGSetRadius(ksp_.implementation(), this->current_radius());
    #else
                KSPCGSetRadius(ksp_.implementation(), this->current_radius());
    #endif
                ierr = KSPSetTolerances(ksp_.implementation(), TRSubproblem::rtol(), TRSubproblem::atol(), PETSC_DEFAULT,  TRSubproblem::max_it());

                ksp_.rtol(TRSubproblem::rtol());
                ksp_.atol(TRSubproblem::atol());
                ksp_.stol(TRSubproblem::stol());
                ksp_.max_it(TRSubproblem::max_it());
                ksp_.verbose(TRSubproblem::verbose());
            }
            else
            {
                // std::cout<<"heeere---- \n";
                PetscErrorCode ierr; UTOPIA_UNUSED(ierr);

                ierr = KSPSetTolerances(ksp_.implementation(), TRSubproblem::rtol(), TRSubproblem::atol(), PETSC_DEFAULT,  TRSubproblem::max_it());

                ksp_.rtol(TRSubproblem::rtol());
                ksp_.atol(TRSubproblem::atol());
                ksp_.stol(TRSubproblem::stol());
                ksp_.max_it(TRSubproblem::max_it());
                ksp_.verbose(TRSubproblem::verbose());
                ksp_.set_initial_guess_non_zero(false);

                PC pc_redundant;
                KSPGetPC(ksp_.implementation(), &pc_redundant);

                KSP innerksp;
                PCRedundantGetKSP(pc_redundant, &innerksp);

                KSPType inner_ksp_type;
                KSPGetType(innerksp, &inner_ksp_type);

                // std::cout<<"inner_ksp_type: "<< inner_ksp_type << "  \n";

                #if UTOPIA_PETSC_VERSION_LESS_THAN(3,8,0)
                    if(std::strcmp(inner_ksp_type,"qcg") == 0)
                        ierr = KSPQCGSetTrustRegionRadius(innerksp, this->current_radius());
                    else if(std::strcmp(inner_ksp_type, "gltr") == 0)
                        ierr = KSPGLTRSetRadius(innerksp, this->current_radius());
                    else if(std::strcmp(inner_ksp_type, "nash") == 0)
                        ierr = KSPNASHSetRadius(innerksp, this->current_radius());
                    else
                        ierr = KSPSTCGSetRadius(innerksp, this->current_radius());
                #else
                    KSPCGSetRadius(innerksp, this->current_radius());
                #endif

            }

        }


    protected:
        KSPSolver ksp_;
        bool redundant_solve_flg_;
    };


    template<typename Matrix, typename Vector>
    class SteihaugToint<Matrix, Vector, PETSC> final: public KSP_TR<Matrix, Vector, PETSC> {

        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "utopia::KSP_TR:: only works with petsc types");

    public:
        SteihaugToint(const std::string &preconditioner = "jacobi")
        : KSP_TR<Matrix, Vector, PETSC>()
        {
            this->ksp_.pc_type(preconditioner);
            this->ksp_.ksp_type("stcg");
        }

        SteihaugToint<Matrix, Vector, PETSC>* clone() const override {
            return new SteihaugToint<Matrix, Vector, PETSC>(*this);
        }

    };

    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class Nash {};

    template<typename Matrix, typename Vector>
    class Nash<Matrix, Vector, PETSC> final: public KSP_TR<Matrix, Vector, PETSC> {

        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "utopia::KSP_TR:: only works with petsc types");

    public:
        Nash(const std::string &preconditioner = "jacobi")
        : KSP_TR<Matrix, Vector, PETSC>()
        {
            this->ksp_.pc_type(preconditioner);
            this->ksp_.ksp_type("nash");
        }

        Nash<Matrix, Vector, PETSC>* clone() const override {
            return new Nash<Matrix, Vector, PETSC>(*this);
        }
    };

    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class Lanczos {};

       template<typename Matrix, typename Vector>
    class Lanczos<Matrix, Vector, PETSC> final: public KSP_TR<Matrix, Vector, PETSC> {

        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "utopia::KSP_TR:: only works with petsc types");

    public:
        Lanczos(const std::string &preconditioner = "jacobi")
        : KSP_TR<Matrix, Vector, PETSC>()
        {
            this->ksp_.pc_type(preconditioner);
            this->ksp_.ksp_type("gltr");
        }

        Lanczos<Matrix, Vector, PETSC>* clone() const override {
            return new Lanczos<Matrix, Vector, PETSC>(*this);
        }

    };

    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class CGNE {};

    template<typename Matrix, typename Vector>
    class CGNE<Matrix, Vector, PETSC> final: public KSP_TR<Matrix, Vector, PETSC> {

        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "utopia::KSP_TR:: only works with petsc types");

    public:
        CGNE(const std::string &preconditioner = "jacobi"): KSP_TR<Matrix, Vector, PETSC>()
        {
            this->ksp_.pc_type(preconditioner);
            this->ksp_.ksp_type("cgne");
        }

        CGNE<Matrix, Vector, PETSC>* clone() const override {
            return new CGNE<Matrix, Vector, PETSC>(*this);
        }

    };

}

#endif //UTOPIA_TR_SUBPROBLEM_KSP_TR_HPP

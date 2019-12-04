#ifndef UTOPIA_BELOS_IMPL_HPP
#define UTOPIA_BELOS_IMPL_HPP

#include "utopia_Belos_solver.hpp"

#include "utopia_make_unique.hpp"
#include "utopia_Wrapper.hpp"

#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>

//TODO remove from here
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Tpetra_CrsMatrix.hpp>

#ifdef WITH_TRILINOS_MUELU
#include <MueLu.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
#else
#warning "Trilinos was not built with MueLu support. AMG cannot be used as a preconditioner for the Belos solver."
#endif //WITH_TRILINOS_MUELU


#ifdef WITH_TRILINOS_IFPACK2
#include <Ifpack2_Factory.hpp>
#else
#warning "Trilinos was not built with Ifpack2 support. Direct preconditioners cannot be used with the Belos solver."
#endif //WITH_TRILINOS_IFPACK2


namespace utopia
{

template <typename Matrix, typename Vector>
class BelosSolver<Matrix, Vector, TRILINOS>::Impl
{
public:
    typedef double ST;

    typedef Tpetra::Operator<>::scalar_type SC;
    typedef Tpetra::Operator<SC>::local_ordinal_type LO;
    typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;

    typedef Kokkos::Compat::KokkosSerialWrapperNode serial_node;

#ifdef  KOKKOS_ENABLE_CUDA
    typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
    typedef cuda_node NT;
#elif defined KOKKOS_ENABLE_ROCM //Kokkos::Compat::KokkosROCmWrapperNode doesn't exist
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::ROCm> rocm_node;
    typedef rocm_node NT;
#elif defined   KOKKOS_ENABLE_OPENMP
    typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
    typedef openmp_node NT;
#else
    typedef serial_node NT;
#endif

    typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
    typedef Tpetra::Operator<SC, LO, GO, NT> OP;

    typedef Belos::LinearProblem<SC, MV, OP> problem_type;
    typedef Belos::SolverManager<SC, MV, OP> solver_type;

#ifdef WITH_TRILINOS_IFPACK2
    typedef Ifpack2::Preconditioner<SC, LO, GO, NT> ifpack_prec_type;
#endif //WITH_TRILINOS_IFPACK2

#ifdef WITH_TRILINOS_MUELU
    typedef MueLu::TpetraOperator<SC, LO, GO, NT> muelu_prec_type;
#endif

    typedef Tpetra::Vector<SC, LO, GO, NT> vec_type;
    typedef Tpetra::CrsMatrix<SC, LO, GO, NT> matrix_type;

    Teuchos::RCP<problem_type> linear_problem;
    Teuchos::RCP<Teuchos::ParameterList> param_list;
    //  auto& utopiaPL;// impl_->param_list->sublist("UTOPIA", true);
    Teuchos::RCP<solver_type> belos_solver;
    Belos::SolverFactory<SC, MV, OP> belos_factory;

    //preconditioner
#ifdef WITH_TRILINOS_IFPACK2
    Teuchos::RCP<ifpack_prec_type> M_ifpack;
#endif //WITH_TRILINOS_IFPACK2

#ifdef WITH_TRILINOS_MUELU
    Teuchos::RCP<muelu_prec_type> M_muelu;
#endif //WITH_TRILINOS_MUELU

};


template <typename Matrix, typename Vector>
BelosSolver<Matrix, Vector, TRILINOS>::BelosSolver(const BelosSolver &other)
    : impl_(utopia::make_unique<Impl>(*other.impl_))
{
    //FIXME
}

template <typename Matrix, typename Vector>
BelosSolver<Matrix, Vector, TRILINOS>::~BelosSolver() {}

template <typename Matrix, typename Vector>
BelosSolver<Matrix, Vector, TRILINOS>::BelosSolver() : impl_(utopia::make_unique<Impl>()) {}

template <typename Matrix, typename Vector>
void BelosSolver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op,
        const std::shared_ptr<const Matrix> &prec)
{
    PreconditionedSolver::update(op, prec);
    // set_problem(*op);
}

template <typename Matrix, typename Vector>
void BelosSolver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op)
{
    PreconditionedSolver::update(op);
    // set_problem(*op);
}

template <typename Matrix, typename Vector>
bool BelosSolver<Matrix, Vector, TRILINOS>::apply(const Vector &rhs, Vector &lhs)
{

    impl_->linear_problem = Teuchos::rcp(
                                new typename Impl::problem_type(
                                    raw_type(*this->get_operator()),
                                    raw_type(lhs),
                                    raw_type(rhs)
                                )
                            );
    set_problem();
    assert(!(impl_->belos_solver.is_null()));
    impl_->belos_solver->solve();
    return true;
}

template <typename Matrix, typename Vector>
int BelosSolver<Matrix, Vector, TRILINOS>::get_num_iter() const
{
    return impl_->belos_solver->getNumIters();
}


template <typename Matrix, typename Vector>
double BelosSolver<Matrix, Vector, TRILINOS>::achieved_tol() const
{
    return impl_->belos_solver->achievedTol();
}

template <typename Matrix, typename Vector>
void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
{
    set_preconditioner(); //(A) //FIXME
}

template <typename Matrix, typename Vector>
void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner(const Matrix &precond)
{
    bool direct_solver = impl_->param_list->template get<bool>("Direct Preconditioner", false);
    std::string dir_prec_type = impl_->param_list->get("Ifpack2 Preconditioner", "prec_type_unset");
    if ( direct_solver )
    {
#ifdef WITH_TRILINOS_IFPACK2
        impl_->M_ifpack = Ifpack2::Factory::create<typename Impl::matrix_type>(dir_prec_type, raw_type(precond));
        assert(!impl_->M_ifpack.is_null());
        impl_->M_ifpack->setParameters(impl_->param_list->sublist(dir_prec_type, false));
        impl_->M_ifpack->initialize();
        impl_->M_ifpack->compute();
        std::string preconditioner_type = impl_->param_list->get("Preconditioner Side", "right");
        std::transform(preconditioner_type.begin(), preconditioner_type.end(), preconditioner_type.begin(), [](unsigned char c)
        {
            return std::tolower(c);
        });
        if (preconditioner_type == "left")
        {
            impl_->linear_problem->setLeftPrec(impl_->M_ifpack);
        }
        else
        {
            impl_->linear_problem->setRightPrec(impl_->M_ifpack);
        }
#else  // WITH_TRILINOS_IFPACK2
        std::cerr << "Cannot use a Direct Preconditioner with the BelosSolver, since Trilinos was not built with Ifpack2 support!" << std::endl;
#endif // WITH_TRILINOS_IFPACK2
    }
    else
    {
#ifdef WITH_TRILINOS_MUELU
        // Multigrid Hierarchy
        impl_->M_muelu = MueLu::CreateTpetraPreconditioner((
                             Teuchos::RCP<typename Impl::OP>) raw_type(precond),
                         impl_->param_list->sublist("MueLu", false)
                                                          );

        assert(!impl_->M_muelu.is_null());
        std::string preconditioner_type = impl_->param_list->get("Preconditioner Side", "right");
        std::transform(preconditioner_type.begin(), preconditioner_type.end(), preconditioner_type.begin(), [](unsigned char c)
        {
            return std::tolower(c);
        });
        if (preconditioner_type == "left")
        {
            impl_->linear_problem->setLeftPrec(impl_->M_muelu);
        }
        else
        {
            impl_->linear_problem->setRightPrec(impl_->M_muelu);
        }
#else
        std::cerr << "Cannot use MueLu as preconditioner since Trilinos was not built with MueLu support." << std::endl;
#endif //WITH_TRILINOS_MUELU
    }
}

template <typename Matrix, typename Vector>
void BelosSolver<Matrix, Vector, TRILINOS>::read_xml(const std::string &path)
{
    if(!path.empty())
    {
        try
        {
            impl_->param_list = Teuchos::getParametersFromXmlFile(path);
        }
        catch(const std::exception &ex)
        {
            std::cerr << ex.what() << std::endl;
            assert(false);
            abort();
        }
    }
    else
    {
        //use default paramlist
    }
}

//read an utopia file and convert in a trilinos list
template <typename Matrix, typename Vector>
void BelosSolver<Matrix, Vector, TRILINOS>::read(Input &in)
{
    Smoother<Matrix, Vector>::read(in);
    PreconditionedSolver::read(in);
    bool enable_adv_opt = false;
    bool enable_flexible_gmres = false;
    std::string solver_type = "CG";

    //in.get("sweeps", sweeps);
    //in.get("relaxation_parameter", _relaxation_parameter);

    //in.get("atol", atol_);
    //in.get("rtol", rtol_);
    //in.get("stol", stol_);
    //in.get("max-it",v max_it_);
    //in.get("verbose", verbose_ );
    impl_->param_list = Teuchos::rcp(new Teuchos::ParameterList);
    in.get("solver_type", solver_type);
    in.get("enable_adv_opt", enable_adv_opt);
    in.get("enable_flexible_gmres",  enable_flexible_gmres );



    in.get("preconditioner", [this](Input &in) {
            std::string precond_type = "prec_type_unset";
            std::string precond_side = "right";
            bool enable_direct_precond = false;
                in.get("enable_direct_precond", enable_direct_precond);
                this->impl_->param_list->set("Direct Preconditioner", enable_direct_precond);
                in.get("precond_side", precond_side);
                in.get("precond_type", precond_type);
                if (enable_direct_precond) {
                    this->impl_->param_list->set("Ifpack2 Preconditioner", precond_type);
                }
                this->impl_->param_list->set("Muelu", precond_type); //TODO create a sublist for Muelu
                this->impl_->param_list->get("Preconditioner Side", precond_side);
              /*this->set_prec(enable_direct_precond);
                this->set_prec(precond_side);
                this->set_prec(precond_type);*/
                });

    impl_->param_list->set("Solver Type", solver_type );
    impl_->param_list->set("Convergence Tolerance", this->rtol());
    impl_->param_list->set("Maximum iterations", this->max_it());
    impl_->param_list->set( "Flexible Gmres", enable_flexible_gmres );

    if (enable_adv_opt)
    {
        int frequency = 10;
        double pol_tol = 1.0e-5;
        int maxrestarts = 20;
        int blocksize = 0;
        int num_blocks = 0;
        int maxdegree = 2;
        bool userandomrhs = false;
        std::string ortho = "ICGS"; /* IMGS is Iterated Modified Gram Schmidt
                                       ICGS is Iterated Classical Gram Schmidt, other is DKGS */
        int recycle = 0;
        bool maxresnorm = true;
        bool estimate_cond_num = true;
        bool assertPositiveDefiniteness = false;
        bool keepHessenberg = false;
        int maxSave = 0;
        int maxDeflate = 0;
        bool use_single_red = false;
        double impTolScale = 0.;
        double impResScale = 0.;
        double expResScale = 0.;
        double expResTest = 0.;
        double damp = 0.;
        double addRoots = 0.;
        double lambda = 0.;
        bool combineConvInner = false;
        int maxiters = 100;

        in.get("frequency", frequency );
        in.get("pol_tol",  pol_tol );
        in.get("maxrestarts",  maxrestarts );
        in.get("blocksize",  blocksize );
        in.get("num_blocks",  num_blocks );
        in.get("maxdegree",  maxdegree );
        in.get("userandomrhs",  userandomrhs );
        in.get("estimate_cond_numortho",  ortho );
        in.get("recycle",  recycle );
        in.get("maxresnorm",  maxresnorm );
        in.get("ortho",  ortho );
        in.get("maxSave",  maxSave );
        in.get("maxDeflate",  maxDeflate );
        in.get("estimate_cond_num",  estimate_cond_num );
        in.get("assertPositiveDefiniteness", assertPositiveDefiniteness  );
        in.get("keepHessenberg",  keepHessenberg );
        in.get("use_single_red",  use_single_red );
        in.get("impTolScale",  impTolScale );
        in.get("impResScale",  impResScale );
        in.get("expResScale",  expResScale );
        in.get("expResTest",  expResTest );
        in.get("damp",  damp );
        in.get("addRoots",  addRoots );
        in.get("lambda",  lambda );
        in.get("combineConvInner", combineConvInner);
        in.get("maxiters", maxiters);

        impl_->param_list->set("Output Frequency", frequency);
        //impl_->param_list->set("S tolerance", this->stol(), "CG");
        //impl_->param_list->set("A tolerance", this->atol(), "CG");
        impl_->param_list->set("Polynomial Tolerance", pol_tol ); // Polynomial convergence tolerance
        impl_->param_list->set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
        impl_->param_list->set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
        impl_->param_list->set( "Num Blocks", num_blocks);             // Maximum number of blocks in Krylov factorization
        impl_->param_list->set( "Orthogonalization", ortho );
        impl_->param_list->set( "Num Recycled Blocks", recycle );
        impl_->param_list->set( "Show Maximum Residual Norm Only", maxresnorm );
        //PCPG
        impl_->param_list->set( "Num Saved Blocks", maxSave );
        impl_->param_list->set( "Num Deflated Blocks", maxDeflate );    // Number of vectors in seed space
        //block GMRES
        //impl_->param_list->set( "Outer Solver", outersolver );
        //impl_->param_list->set( "Outer Solver Params", belosList ); //two nested belos solvers
        impl_->param_list->set( "Estimate Condition Number", estimate_cond_num );
        impl_->param_list->set( "Use Single Reduction", use_single_red ); // Use single reduction CG iteration
        impl_->param_list->set( "Fold Convergence Detection Into Allreduce", combineConvInner );
        impl_->param_list->set("Keep Hessenberg", keepHessenberg);
        impl_->param_list->set("Implicit Tolerance Scale Factor", impTolScale );
        impl_->param_list->set("Implicit Residual Scaling", impResScale );
        impl_->param_list->set("Explicit Residual Scaling", expResScale );
        impl_->param_list->set("Explicit Residual Test", expResTest);
        //impl_->param_list->set("Deflation Quorum", static_cast<int>(defQuorum_default_)
        //impl_->param_list->set("Adaptive Block Size", static_cast<bool>(adaptiveBlockSize_default_),
        //impl_->param_list->set("Show Maximum Residual Norm Only", static_cast<bool>(showMaxResNormOnly_default_),
        //impl_->param_list->set("Implicit Residual Scaling", resScale_default_,
        impl_->param_list->set("Damped Poly", damp);
        impl_->param_list->set("Add Roots", addRoots);
        impl_->param_list->set("Assert Positive Definiteness",assertPositiveDefiniteness);
        impl_->param_list->set("Max Size For Condest",maxiters);
        impl_->param_list->set( "Lambda", lambda );
        impl_->param_list->set( "Maximum Degree", maxdegree );          // Maximum degree of the GMRES polynomial
        impl_->param_list->set( "Random RHS", userandomrhs );
    }

    if (this->verbose())
    {
        impl_->param_list->set("Verbose", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails + Belos::FinalSummary);
        impl_->param_list->set("Output Style", Belos::General);
    }
    else
    {
        impl_->param_list->set("Verbose", Belos::Errors);
        impl_->param_list->set("Output Style", Belos::Brief);
    }

    if (this->verbose())
    {
        std::cout << "Current Parameters:" << std::endl;
        impl_->param_list->print();
        if (impl_->belos_solver.is_null())
        {std::cout << "Belos solver is currently null" << std::endl;}
    }


}

// available parameters
template <typename Matrix, typename Vector>
void BelosSolver<Matrix, Vector, TRILINOS>::print_usage(std::ostream &os ) const
{
    Smoother<Matrix, Vector>::print_usage(os);
    PreconditionedSolver::print_usage(os);
    impl_->param_list->print();
    std::cout << "Valid solver Parameters:" << std::endl;
    impl_->belos_solver->getValidParameters()->print();
    std::cout << "Current solver Parameters:" << std::endl;
    impl_->belos_solver->getCurrentParameters()->print();
}

template <typename Matrix, typename Vector>
BelosSolver<Matrix, Vector, TRILINOS> * BelosSolver<Matrix, Vector, TRILINOS>::clone() const
{
    return new BelosSolver(*this);
}

template <typename Matrix, typename Vector>
bool BelosSolver<Matrix, Vector, TRILINOS>::smooth(const Vector &rhs, Vector &x)
{
    return false;
}

template <typename Matrix, typename Vector>
bool BelosSolver<Matrix, Vector, TRILINOS>::set_problem()
{
    impl_->linear_problem->setProblem();
    auto sol_type = impl_->param_list->get("Solver Type", "CG");
    auto belos_params = Teuchos::sublist(impl_->param_list, sol_type, false);
    impl_->belos_solver = impl_->belos_factory.create( sol_type, belos_params); //to change it to have the specialization
    impl_->belos_solver->setProblem(impl_->linear_problem);
    if (this->verbose())
    {
        std::cout << "Current Parameters when setting the Problem" << std::endl;
        impl_->belos_solver->getCurrentParameters()->print();
    }
    set_preconditioner();
    return true;
}

template <typename Matrix, typename Vector>
bool BelosSolver<Matrix, Vector, TRILINOS>::set_problem(Matrix &A)
{
    impl_->linear_problem->setProblem();
    auto sol_type = impl_->param_list->get("Solver Type", "CG");
    auto belos_params = Teuchos::sublist(impl_->param_list, sol_type, false);
    impl_->belos_solver = impl_->belos_factory.create( sol_type, belos_params); //to change it to have the specialization
    set_preconditioner(); //(A);
    impl_->belos_solver->setProblem(impl_->linear_problem);
    if (this->verbose())
    {
        std::cout << std::endl << "Actual Parameters used by Belos:" << std::endl;
        impl_->belos_solver->getCurrentParameters()->print();
    }
    return true;
}

template <typename Matrix, typename Vector>
void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner()//const std::shared_ptr<Preconditioner> &precond)
{
    bool direct_solver = impl_->param_list->template get<bool>("Direct Preconditioner", false);
    std::string dir_prec_type = impl_->param_list->get("Ifpack2 Preconditioner", "prec_type_unset");

    if ( direct_solver )
    {
#ifdef WITH_TRILINOS_IFPACK2
        impl_->M_ifpack = Ifpack2::Factory::create<typename Impl::matrix_type>(dir_prec_type, raw_type(*this->get_operator()));
        assert(!impl_->M_ifpack.is_null());
        impl_->M_ifpack->setParameters(impl_->param_list->sublist(dir_prec_type, false));
        impl_->M_ifpack->initialize();
        impl_->M_ifpack->compute();
        std::string preconditioner_side = impl_->param_list->get("Preconditioner Side", "right");
        std::transform(preconditioner_side.begin(), preconditioner_side.end(), preconditioner_side.begin(), [](unsigned char c)
        {
            return std::tolower(c);
        });
        if (preconditioner_side == "left")
        {
            impl_->linear_problem->setLeftPrec(impl_->M_ifpack);
        }
        else
        {
            impl_->linear_problem->setRightPrec(impl_->M_ifpack);
        }
#else  //WITH_TRILINOS_IFPACK2
        std::cerr << "Cannot use a Direct Preconditioner with the BelosSolver, since Trilinos was not built with Ifpack2 support!" << std::endl;
#endif //WITH_TRILINOS_IFPACK2
    }
    else
    {
#ifdef WITH_TRILINOS_MUELU
        // Multigrid Hierarchy
        impl_->M_muelu = MueLu::CreateTpetraPreconditioner(raw_type(*this->get_operator()), impl_->param_list->sublist("MueLu", false));
        assert(!impl_->M_muelu.is_null());
        std::string preconditioner_side = impl_->param_list->get("Preconditioner Side", "right");
        std::transform(preconditioner_side.begin(), preconditioner_side.end(), preconditioner_side.begin(), [](unsigned char c)
        {
            return std::tolower(c);
        });
        if (preconditioner_side == "left")
        {
            impl_->linear_problem->setLeftPrec(impl_->M_muelu);
        }
        else
        {
            impl_->linear_problem->setRightPrec(impl_->M_muelu);
        }
#else  // WITH_TRILINOS_MUELU
        std::cerr << "Cannot use MueLu as preconditioner since Trilinos was not built with MueLu support." << std::endl;
#endif //WITH_TRILINOS_MUELU
    }
}

}  // namespace utopia

#endif //UTOPIA_BELOS_IMPL_HPP

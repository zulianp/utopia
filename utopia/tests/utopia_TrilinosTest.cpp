#include "utopia_Base.hpp"
//#include "utopia_trilinos_Vector.hpp"
#ifdef WITH_TRILINOS
#include "utopia_TrilinosTest.hpp"
#include "utopia.hpp"
#include "utopia_trilinos_solvers.hpp"
#include <iostream>

using namespace Teuchos;

void trilinos_example_test()
    {
    typedef Tpetra::Map<>                         map_type;
    typedef Tpetra::Vector<>::scalar_type         scalar_type;
    typedef Tpetra::Vector<>::local_ordinal_type  local_ordinal_type;
    typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
    typedef Tpetra::Vector<>::mag_type            magnitude_type;
    typedef Tpetra::CrsMatrix<>                   crs_matrix_type;

    // The number of rows and columns in the matrix.
    const Tpetra::global_size_t numGblIndices = 50;
    const global_ordinal_type indexBase = 0;

    // Number of iterations
    const int niters = 500;
    // Desired (absolute) residual tolerance
    const magnitude_type tolerance = 1.0e-2;


    Teuchos::RCP<const Comm<int> > comm = rcp(new MpiComm<int>(MPI_COMM_WORLD));

    // Make an empty new parameter list.
    Teuchos::RCP<Teuchos::ParameterList> solverParams = parameterList();

    if (comm->getRank () == 0)
        {
        // On (MPI) Process 0, print out the Tpetra software version.
        std::cout << Tpetra::version () << std::endl << std::endl;
        }

    // create the Map
    Teuchos::RCP<const map_type> map = rcp (new map_type (numGblIndices, indexBase, comm));
    const size_t numMyElements = map->getNodeNumElements ();


    std::cout<< "Creating the sparse matrix **" << std::endl;
    // Create a Tpetra sparse matrix whose rows have distribution given by the Map.
    //Teuchos::RCP<crs_matrix_type> A (new crs_matrix_type (map, indexBase));

    std::cout<< "1" << std::endl;
    // Fill the sparse matrix, one row at a time.
    const scalar_type two    = static_cast<scalar_type> (2.0);
    const scalar_type negOne = static_cast<scalar_type> (-1.0);
    Teuchos::RCP<const Teuchos::Comm<int> > comm2 = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));




    /*    DSMatrixd A;
        assemble_laplacian_1D(n, A);

        // exact solution to our problem
        const DVectord u_exact  = values(n, solution);
        //
        // constructing initial guess
        DVectord u = zeros(n);
        //
        // constructing rhs
        const DVectord rhs   = A * u_exact;
        //
        // setting up parameters of solver
        Parameters params;
        params.tol(1e-9);
        params.lin_solver_type("UTOPIA_CG");
        params.linear_solver_verbose(true);
        //
        // solve
        solve(A, rhs, u, params);
    */
    std::cout<< "2" << std::endl;
    utopia::TpetraMatrix A_void2();
    utopia::TpetraMatrix A_void(comm2);
    utopia::TpetraMatrix A_map(map);
    utopia::TpetraMatrix A_mat(A_map);
    std::cout<< "3" << std::endl;
    utopia::TpetraVector x_void();
    utopia::TpetraVector x_map(map);
    std::cout<< "4" << std::endl;
//size_t numMyElements_void = A_void2.getNodeNumElements();

    for (local_ordinal_type lclRow = 0; lclRow < static_cast<local_ordinal_type> (numMyElements); ++lclRow)
        {
        const global_ordinal_type gblRow = map->getGlobalElement (lclRow);
        //A->insertGlobalValues (gblRow, Teuchos::tuple<global_ordinal_type> (gblRow, gblRow + 1), Teuchos::tuple<scalar_type> (two, negOne));
        A_map.insertGlobalValues (gblRow, Teuchos::tuple<global_ordinal_type> (gblRow, gblRow + 1), Teuchos::tuple<scalar_type> (two, negOne));
        }

    std::cout<< "5" << std::endl;

    const size_t numMyElements_void = A_void.getNodeNumElements();

    std::cout<< "1a" << std::endl;
    for (local_ordinal_type lclRow = 0; lclRow < static_cast<local_ordinal_type> (numMyElements_void); ++lclRow)
        {
        std::cout<< "a" << lclRow << std::endl;
        const global_ordinal_type gblRow = A_void.getGlobalElement (lclRow);
        std::cout<< "b" << std::endl;
        A_void.insertGlobalValues (gblRow, Teuchos::tuple<global_ordinal_type> (gblRow, gblRow + 1), Teuchos::tuple<scalar_type> (two, negOne));
        std::cout<< "c" << std::endl;
        }
    std::cout<< "6" << std::endl;
    // Tell the sparse matrix that we are done adding entries to it.
    //A->fillComplete();
    std::cout<< "7" << std::endl;
    A_void.fillComplete();
    std::cout<< "8" << std::endl;
    A_map.fillComplete();
    std::cout<< "9" << std::endl;

    // Run the power method and report the result.
    scalar_type lambda;// = PowerMethod<crs_matrix_type>::run (*A_map, niters, tolerance, out);

    std::cout << "Estimated max eigenvalue: " << lambda << std::endl;
//////////////////////////////////////////////////////////////////////////////////////




    typedef double                           ST;
    typedef Teuchos::ScalarTraits<ST>                SCT;
    typedef SCT::magnitudeType               MT;
    typedef Tpetra::Operator<ST,int>         OP;
    typedef Tpetra::MultiVector<ST,int>      MV;
    typedef Belos::OperatorTraits<ST,MV,OP> OPT;
    typedef Belos::MultiVecTraits<ST,MV>    MVT;

//    Teuchos::GlobalMPISession mpisess(&argc,&argv,&cout);

    bool success = false;
    bool verbose = false;
   
        {
        const ST one  = SCT::one();

        int MyPID = 0;

        typedef Tpetra::DefaultPlatform::DefaultPlatformType           Platform;
        typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

        Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
        RCP<const Comm<int> > comm = platform.getComm();
        RCP<Node>             node = platform.getNode();

        //
        // Get test parameters from command-line processor
        //
        bool proc_verbose = false;
        bool debug = false;
        int frequency = -1;  // how often residuals are printed by solver
        int numrhs = 1;      // total number of right-hand sides to solve for
        int blocksize = 1;   // blocksize used by solver
        int maxiters = -1;   // maximum number of iterations for solver to use
        std::string filename("bcsstk14.hb");
        MT tol = 1.0e-5;     // relative residual tolerance

        CommandLineProcessor cmdp(false,true);

        cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
        cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
        cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
        cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
        cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
        cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
        cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
        cmdp.setOption("block-size",&blocksize,"Block size to be used by the CG solver.");
        /*if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL)
            {
            return -1;
            }
        if (debug)
            {
            verbose = true;
            }
        if (!verbose)
            {
            frequency = -1;  // reset frequency if test is not verbose
            }
        */
        MyPID = rank(*comm);
        proc_verbose = ( verbose && (MyPID==0) );

        if (proc_verbose)
            {
            std::cout << Belos::Belos_Version() << std::endl << std::endl;
            }

        //
        // Get the data from the HB file and build the Map,Matrix
        //
        Teuchos::RCP<Tpetra::CrsMatrix<ST,int> > A;
//        Tpetra::Utils::readHBMatrix(filename,comm,node,A);
        Teuchos::RCP<const Tpetra::Map<int> > map = A->getDomainMap();

        // Create initial vectors
        Teuchos::RCP<Tpetra::MultiVector<ST,int> > B, X;
        X = rcp( new Tpetra::MultiVector<ST,int>(map,numrhs) );
        MVT::MvRandom( *X );
        B = rcp( new Tpetra::MultiVector<ST,int>(map,numrhs) );


        OPT::Apply( *A, *X, *B );
        MVT::MvInit( *X, 0.0 );

        //
        // ********Other information used by block solver***********
        // *****************(can be user specified)******************
        //
        const int NumGlobalElements = B->getGlobalLength();
        if (maxiters == -1)
            {
            maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
            }
        //
        ParameterList belosList;
        belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
        belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
        belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
        int verbLevel = Belos::Errors + Belos::Warnings;
        if (debug)
            {
            verbLevel += Belos::Debug;
            }
        if (verbose)
            {
            verbLevel += Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
            }
        belosList.set( "Verbosity", verbLevel );
        if (verbose)
            {
            if (frequency > 0)
                {
                belosList.set( "Output Frequency", frequency );
                }
            }
        //
        // Construct an unpreconditioned linear problem instance.
        //
        Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
        bool set = problem.setProblem();
        if (set == false)
            {
            if (proc_verbose)
                std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
            }


// *************Start the block CG iteration***********************
        // *******************************************************************
        //
        Belos::BlockCGSolMgr<ST,MV,OP> solver( rcp(&problem,false), rcp(&belosList,false) );

        //
        // **********Print out information about problem*******************
        //
        if (proc_verbose)
            {
            std::cout << std::endl << std::endl;
            std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
            std::cout << "Number of right-hand sides: " << numrhs << std::endl;
            std::cout << "Block size used by solver: " << blocksize << std::endl;
            std::cout << "Max number of CG iterations: " << maxiters << std::endl;
            std::cout << "Relative residual tolerance: " << tol << std::endl;
            std::cout << std::endl;
            }
        //
        // Perform solve
        //
        Belos::ReturnType ret = solver.solve();
        //
        // Compute actual residuals.
        //
        bool badRes = false;
        std::vector<MT> actual_resids( numrhs );
        std::vector<MT> rhs_norm( numrhs );
        Tpetra::MultiVector<ST,int> resid(map, numrhs);
        OPT::Apply( *A, *X, resid );
        MVT::MvAddMv( -one, resid, one, *B, resid );
        MVT::MvNorm( resid, actual_resids );
        MVT::MvNorm( *B, rhs_norm );
        if (proc_verbose)
            {

            std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
            }
        for ( int i=0; i<numrhs; i++)
            {
            MT actRes = actual_resids[i]/rhs_norm[i];
            if (proc_verbose)
                {
                std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
                }
            if (actRes > tol) badRes = true;
            }

        success = (ret==Belos::Converged && !badRes);

    /*    if (success)
            {
            if (proc_verbose)
                std::cout << "\nEnd Result: TEST PASSED" << std::endl;
            }
        else
            {
            if (proc_verbose)
                std::cout << "\nEnd Result: TEST FAILED" << std::endl;
            }*/
//


//
        }
    }
}
namespace utopia
{
void run_trilinos_test()
    {
    UTOPIA_UNIT_TEST_BEGIN("TrilinosTest");
    UTOPIA_RUN_TEST(trilinos_example_test);
    UTOPIA_UNIT_TEST_END("TrilinosTest");
    }
}

#else //WITH_TRILINOS
namespace utopia
{
void run_trilinos_test() {}



// A(0, 0:1) = [2, -1]
/*  if (gblRow == 0) { A->insertGlobalValues (gblRow, Teuchos::tuple<global_ordinal_type> (gblRow, gblRow + 1), Teuchos::tuple<scalar_type> (two, negOne));
  }
  // A(N-1, N-2:N-1) = [-1, 2]
  else if (static_cast<Tpetra::global_size_t> (gblRow) == numGblIndices - 1)
    {A->insertGlobalValues (gblRow,Teuchos::tuple<global_ordinal_type> (gblRow - 1, gblRow),Teuchos::tuple<scalar_type> (negOne, two));
  }
  // A(i, i-1:i+1) = [-1, 2, -1]
  else {A->insertGlobalValues (gblRow,Teuchos::tuple<global_ordinal_type> (gblRow - 1, gblRow, gblRow + 1),Teuchos::tuple<scalar_type> (negOne, two, negOne));
  }*/
}
#endif //WITH_TRILINOS

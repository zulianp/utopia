

//
// Example: Create an Ifpack2 preconditioner from a Tpetra::CrsMatrix.
//

#include <Ifpack2_Factory.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <iostream>

// The PreconditionerFactory class encapsulates creation of an Ifpack2
// preconditioner (returned as a Tpetra::Operator, rather than an
// Ifpack2::Preconditioner) from a Tpetra::CrsMatrix.  Change the
// preconditioner parameters by changing the function that creates a
// ParameterList for Ifpack2.  See the Ifpack2 documentation for a
// list and description of the parameters that it accepts.
//
// We template PreconditionerFactory on the matrix type, rather than
// on the template arguments of Tpetra::CrsMatrix, because CrsMatrix
// has some nice typedefs that let us retrieve those template
// arguments.  It's easier to template on one thing than on five
// things!
template<class TpetraMatrixType>
class PreconditionerFactory
    {
    public:
        // Fetch the typedefs defined by Tpetra::CrsMatrix.
        typedef typename TpetraMatrixType::scalar_type scalar_type;
        typedef typename TpetraMatrixType::local_ordinal_type local_ordinal_type;
        typedef typename TpetraMatrixType::global_ordinal_type global_ordinal_type;
        typedef typename TpetraMatrixType::node_type node_type;

        // Making this typedef public lets users of PreconditionerFactory
        // get straight to the original template argument.
        typedef TpetraMatrixType matrix_type;

        // A Tpetra::Operator is an abstraction of a function mapping a
        // (Multi)Vector to a (Multi)Vector.
        typedef Tpetra::Operator<scalar_type, local_ordinal_type,
                global_ordinal_type, node_type> op_type;

    private:
        // These are just some convenience typedefs.
        typedef Teuchos::ScalarTraits<scalar_type> STS;
        typedef typename STS::magnitudeType magnitude_type;
        typedef Teuchos::ScalarTraits<magnitude_type> STM;

        // An Ifpack2::Preconditioner is-a Tpetra::Operator.  Ifpack2
        // creates a Preconditioner object, but users of iterative methods
        // want a Tpetra::Operator.  That's why create() returns an
        // Operator instead of a Preconditioner.
        typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type,
                global_ordinal_type, node_type> prec_type;

    public:

        // The constructor doesn't do anything, since this factory doesn't
        // keep any state.
        PreconditionerFactory () {}

        // Return a ParameterList for asking Ifpack2 to create an ILUT
        // incomplete factorization preconditioner with fill level 2, drop
        // tolerance 0, and absolute threshold 0.1.
        Teuchos::RCP<Teuchos::ParameterList>
        parameterListForIfpack2 () const
            {
            using Teuchos::ParameterList;
            using Teuchos::parameterList;
            using Teuchos::RCP;

            // The name of the type of preconditioner to use.
            const std::string precondType ("ILUT");

            // Ifpack2 expects double-precision arguments here.
            const double fillLevel = 2.0;
            const double dropTol = 0.0;
            const double absThreshold = 0.1;
            const bool verbose = true;
            const bool debug = false;

            RCP<ParameterList> pl = parameterList ("Preconditioner");
            pl->set ("Ifpack2::Preconditioner", precondType);

            ParameterList precParams ("Ifpack2");
            precParams.set ("fact: ilut level-of-fill", fillLevel);
            precParams.set ("fact: drop tolerance", dropTol);
            precParams.set ("fact: absolute threshold", absThreshold);

            pl->set ("Ifpack2", precParams);
            return pl;
            }

        // Compute and return an Ifpack2 preconditioner.
        Teuchos::RCP<op_type>
        create (const Teuchos::RCP<const matrix_type>& A,
                const Teuchos::RCP<Teuchos::ParameterList>& plist,
                std::ostream& out,
                std::ostream& err) const
            {
            using Teuchos::outArg;
            using Teuchos::ParameterList;
            using Teuchos::parameterList;
            using Teuchos::RCP;
            using Teuchos::rcp;
            using Teuchos::Time;
            using Teuchos::TimeMonitor;
            using std::endl;

            RCP<Time> initTimer = TimeMonitor::getNewCounter ("Ifpack2::Preconditioner::initialize");
            RCP<Time> computeTimer = TimeMonitor::getNewCounter ("Ifpack2::Preconditioner::compute");
            RCP<Time> condestTimer = TimeMonitor::getNewCounter ("Ifpack2::Preconditioner::condest");

            err << "Creating ILUT preconditioner" << endl
                << "-- Configuring" << endl;
            RCP<prec_type> prec;
                {
                Ifpack2::Factory factory;

                // Get the preconditioner type.
                const std::string precName =
                    plist->get<std::string> ("Ifpack2::Preconditioner");

                // Set up the preconditioner of that type.
                prec = factory.create (precName, A);

                ParameterList ifpack2Params;
                if (plist->isSublist ("Ifpack2"))
                    ifpack2Params = plist->sublist ("Ifpack2");
                else
                    ifpack2Params.setName ("Ifpack2");
                prec->setParameters (ifpack2Params);
                }

            err << "-- Initializing" << endl;
                {
                TimeMonitor mon (*initTimer);
                prec->initialize();
                }
            err << "-- Computing" << endl;
                {
                TimeMonitor mon (*computeTimer);
                prec->compute();
                }

            err << "-- Estimating condition number" << endl;
            magnitude_type condest = STM::one();
                {
                TimeMonitor mon (*condestTimer);
                condest = prec->computeCondEst (Ifpack2::Cheap);
                }

            return prec;
            }
    };

// Create and return a simple example CrsMatrix.
template<class TpetraMatrixType>
Teuchos::RCP<const TpetraMatrixType>
createMatrix (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
              const Teuchos::RCP<typename TpetraMatrixType::node_type>& node)
    {
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    using Teuchos::tuple;

    typedef TpetraMatrixType matrix_type;

    // Fetch the timer for sparse matrix creation.
    //
    // If you are using Trilinos 10.6 instead of the development branch
    // (10.7), just create a new timer here, and remove the bit in
    // main() that creates a timer.
    //
    RCP<Time> timer = TimeMonitor::lookupCounter ("Sparse matrix creation");
    if (timer.is_null())
        timer = TimeMonitor::getNewCounter ("Sparse matrix creation");

    // Time the whole scope of this routine, not counting timer lookup.
    TimeMonitor monitor (*timer);

    // Fetch typedefs from the Tpetra::CrsMatrix.
    typedef typename TpetraMatrixType::scalar_type scalar_type;
    typedef typename TpetraMatrixType::local_ordinal_type local_ordinal_type;
    typedef typename TpetraMatrixType::global_ordinal_type global_ordinal_type;
    typedef typename TpetraMatrixType::node_type node_type;

    // The type of the Tpetra::Map that describes how the matrix is distributed.
    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    // The global number of rows in the matrix A to create.  We scale
    // this relative to the number of (MPI) processes, so that no matter
    // how many MPI processes you run, every process will have 10 rows.
    const Tpetra::global_size_t numGlobalElements = 10 * comm->getSize();

    // Construct a Map that puts approximately the same number of
    // equations on each processor.
    const global_ordinal_type indexBase = 0;
    RCP<const map_type > map =
        rcp (new map_type (numGlobalElements, indexBase, comm,
                           Tpetra::GloballyDistributed, node));

    // Get update list and the number of equations that this MPI process
    // owns.
    const size_t numMyElements = map->getNodeNumElements();
    ArrayView<const global_ordinal_type> myGlobalElements = map->getNodeElementList();

    // NumNz[i] will be the number of nonzero entries for the i-th
    // global equation on this MPI process.
    ArrayRCP<size_t> NumNz = arcp<size_t> (numMyElements);

    // We are building a tridiagonal matrix where each row is (-1 2 -1),
    // so we need 2 off-diagonal terms (except for the first and last
    // equation).
    for (size_t i = 0; i < numMyElements; ++i)
        {
        if (myGlobalElements[i] == 0 || static_cast<Tpetra::global_size_t>(myGlobalElements[i]) == numGlobalElements-1)
            {
            NumNz[i] = 2; // First or last equation
            }
        else
            {
            NumNz[i] = 3;
            }
        }

    // Create a Tpetra::Matrix using the Map, with a static allocation
    // dictated by NumNz.
    RCP<matrix_type> A = rcp (new matrix_type (map, NumNz, Tpetra::StaticProfile));

    // We are done with NumNZ; free it.
    NumNz = Teuchos::null;

    // Add rows one at a time.  Off diagonal values will always be -1.
    const scalar_type two    = static_cast<scalar_type>( 2.0);
    const scalar_type negOne = static_cast<scalar_type>(-1.0);

    for (size_t i = 0; i < numMyElements; i++)
        {
        if (myGlobalElements[i] == 0)
            {
            A->insertGlobalValues (myGlobalElements[i],
                                   tuple (myGlobalElements[i], myGlobalElements[i]+1),
                                   tuple (two, negOne));
            }
        else if (static_cast<Tpetra::global_size_t> (myGlobalElements[i]) == numGlobalElements-1)
            {
            A->insertGlobalValues (myGlobalElements[i],
                                   tuple (myGlobalElements[i]-1, myGlobalElements[i]),
                                   tuple (negOne, two));
            }
        else
            {
            A->insertGlobalValues (myGlobalElements[i],
                                   tuple (myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1),
                                   tuple (negOne, two, negOne));
            }
        }

    // Finish up the matrix.
    A->fillComplete ();
    return A;
    }


template<class NodeType>
void
example (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const Teuchos::RCP<NodeType>& node,
         std::ostream& out,
         std::ostream& err)
    {
    using std::endl;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Print out the Tpetra software version information.
    out << Tpetra::version() << endl << endl;

    // Set up Tpetra typedefs.
    typedef double scalar_type;
    typedef int local_ordinal_type;
    typedef long global_ordinal_type;
    typedef NodeType node_type;
    typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> matrix_type;
    typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> op_type;

    // Create an example sparse matrix.
    RCP<const matrix_type> A = createMatrix<matrix_type> (comm, node);

    // Create a factory for creating a preconditioner.
    PreconditionerFactory<matrix_type> factory;
    // ParameterList for creating an Ifpack2 preconditioner.
    RCP<ParameterList> plist = factory.parameterListForIfpack2 ();
    // Compute the preconditioner using the matrix A.
    // The matrix A itself is not modified.
    RCP<op_type> Prec = factory.create (A, plist, out, err);
    }


int main (int argc, char *argv[])
    {
    using std::endl;
    using Teuchos::RCP;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;

    Teuchos::oblackholestream blackHole;
    Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
    RCP<const Teuchos::Comm<int> > comm =
        Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

    typedef Kokkos::DefaultNode::DefaultNodeType node_type;
    RCP<node_type> node = Kokkos::DefaultNode::getDefaultNode ();

    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();
    std::ostream& out = (myRank == 0) ? std::cout : blackHole;
    std::ostream& err = (myRank == 0) ? std::cerr : blackHole;

    // Make a timer for sparse matrix creation.
    //
    // If you are using Trilinos 10.6 instead of the development branch
    // (10.7), just delete this line of code, and make the other change
    // mentioned above.
    RCP<Time> sparseMatrixCreationTimer =
        TimeMonitor::getNewCounter ("Sparse matrix creation");

    // Run the whole example: create the sparse matrix, and compute the preconditioner.
    example (comm, node, out, err);

    // Summarize global performance timing results, for all timers
    // created using TimeMonitor::getNewCounter().
    TimeMonitor::summarize (out);

    return 0;
    }

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
    return -1;
    }
//
// *******************************************************************
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
MultiVector<ST,int> resid(map, numrhs);
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

if (ret!=Belos::Converged || badRes)
    {
    if (proc_verbose)
        {
        std::cout << "\nEnd Result: TEST FAILED" << std::endl;
        }
    return -1;
    }
//
// Default return value
//
if (proc_verbose)
    {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
    }
return 0;
//
    } // end test_bl_cg_hb.cp

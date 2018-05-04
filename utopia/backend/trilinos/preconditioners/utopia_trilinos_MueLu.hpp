M
UE
L
U
as a preconditioner within B
ELOS
The following code shows the basic steps of how to use a M
UE
L
U
multigrid preconditioner
with T
PETRA
linear algebra library and with a linear solver from B
ELOS
.  To keep the example
short and clear, we skip the template parameters and focus on the algorithmic outline of setting up
a linear solver. For further details, a user may refer to the
examples
and
test
directories.
First, we create the M
UE
L
U
multigrid preconditioner.  It can be done in a few ways.  For in-
stance, multigrid parameters can be read from an XML file (e.g.,
        mueluOptions.xml
        in the example
        below).
Teuchos::RCP<Tpetra::CrsMatrix<> > A;
// create A here ...
std::string optionsFile = "mueluOptions.xml";
Teuchos::RCP<MueLu::TpetraOperator> mueLuPreconditioner =
    MueLu::CreateTpetraPreconditioner(A, optionsFile);
The XML file contains multigrid options.  A typical file with M
UE
L
U
parameters looks like the
following.
<ParameterList name="MueLu">
                    <Parameter name="verbosity" type="string" value="low"/>
                                    <Parameter name="max levels" type="int" value="3"/>
                                            <Parameter name="coarse: max size" type="int" value="10"/>
                                                    <Parameter name="multigrid algorithm" type="string" value="sa"/>
                                                            <!-- Damped Jacobi smoothing -->
                                                            <Parameter name="smoother: type" type="string" value="RELAXATION"/>
                                                                    <ParameterList name="smoother: params">
                                                                            <Parameter name="relaxation: type" type="string" value="Jacobi"/>
                                                                                    <Parameter name="relaxation: sweeps" type="int" value="1"/>
                                                                                            <Parameter name="relaxation: damping factor" type="double" value="0.9"/>
                                                                                                    </ParameterList>
                                                                                                    <!-- Aggregation -->
                                                                                                    19
                                                                                                    <Parameter name="aggregation: type" type="string" value="uncoupled"/>
                                                                                                            <Parameter name="aggregation: min agg size" type="int" value="3"/>
                                                                                                                    <Parameter name="aggregation: max agg size" type="int" value="9"/>
                                                                                                                            </ParameterList>
                                                                                                                            It defines a three level smoothed aggregation multigrid algorithm. The aggregation size is between
                                                                                                                            three and nine(2D)/27(3D) nodes.  One sweep with a damped Jacobi method is used as a level
                                                                                                                            smoother.  By default, a direct solver is applied on the coarsest level.  A complete list of available
                                                                                                                            parameters and valid parameter choices is given in
                                                                                                                            §
                                                                                                                            5 of this User’s Guide.
                                                                                                                            Users can also construct a multigrid preconditioner using a provided
                                                                                                                            ParameterList
                                                                                                                            without
                                                                                                                            accessing any files in the following manner.
                                                                                                                            Teuchos::RCP<Tpetra::CrsMatrix<> > A;
// create A here ...
Teuchos::ParameterList paramList;
paramList.set("verbosity", "low");
paramList.set("max levels", 3);
paramList.set("coarse: max size", 10);
paramList.set("multigrid algorithm", "sa");
// ...
Teuchos::RCP<MueLu::TpetraOperator> mueLuPreconditioner =
    MueLu::CreateTpetraPreconditioner(A, paramList);
Besides the linear operator
A
, we also need an initial guess vector for the solution
X
and a right
hand side vector
B
for solving a linear system.
    Teuchos::RCP<const Tpetra::Map<> > map = A->getDomainMap();
// Create initial vectors
Teuchos::RCP<Tpetra::MultiVector<> > B, X;
X = Teuchos::rcp( new Tpetra::MultiVector<>(map,numrhs) );
Belos::MultiVecTraits<>::MvRandom( *X );
B = Teuchos::rcp( new Tpetra::MultiVector<>(map,numrhs) );
Belos::OperatorTraits<>::Apply( *A, *X, *B );
Belos::MultiVecTraits<>::MvInit( *X, 0.0 );
To generate a dummy example, the above code first declares two vectors.  Then, a right hand side
vector is calculated as the matrix-vector product of a random vector with the operator
A
.  Finally,
an initial guess is initialized with zeros.
Then, one can define a
Belos::LinearProblem
object where the
mueLuPreconditioner
is
used for left preconditioning
Belos::LinearProblem<> problem( A, X, B );
problem->setLeftPrec(mueLuPreconditioner);
20
bool set = problem.setProblem();
Next, we set up a B
ELOS
solver using some basic parameters
Teuchos::ParameterList belosList;
belosList.set( "Block Size", 1 );
belosList.set( "Maximum Iterations", 100 );
belosList.set( "Convergence Tolerance", 1e-10 );
belosList.set( "Output Frequency", 1 );
belosList.set( "Verbosity", Belos::TimingDetails + Belos::FinalSummary );
Belos::BlockCGSolMgr<> solver( rcp(&problem,false), rcp(&belosList,false) );
Finally, we solve the system.
Belos::ReturnType ret = solver.solve();


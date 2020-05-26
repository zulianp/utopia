/*   Example building stiffness matrix and right hand side for a Poisson equation
 using nodal (Hgrad) elements.

 div grad u = f in Omega
 u = 0 on Gamma

 Discrete linear system for nodal coefficients(x):
 Kx = b

 K - HGrad stiffness matrix
 b - right hand side vector

 int NX              - num intervals in x direction (assumed box domain, 0,1)
 int NY              - num intervals in y direction (assumed box domain, 0,1)
 int NZ              - num intervals in z direction (assumed box domain, 0,1)
 verbose (optional)  - any character, indicates verbose output
 */

// Intrepid2 includes
#include <Intrepid2_ArrayTools.hpp>
#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>
#include <Intrepid2_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid2_RealSpaceTools.hpp>
#include <Intrepid2_Utils.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

// Teuchos includes
#include <Teuchos_BLAS.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_oblackholestream.hpp>

// Belos includes
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>

// Ifpack2 includes
#include <Ifpack2_Factory.hpp>

// Amesos2 includes
#include <Amesos2.hpp>

// Muelu includes
#include <MueLu.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_Utilities.hpp>

// Shards includes
#include "Shards_CellTopology.hpp"

// Tpetra includes
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>

// Kokkos and Tpetra typedefs
// typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
// typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
typedef Kokkos::Compat::KokkosSerialWrapperNode serial_node;
// typedef Kokkos::Compat::KokkosThreadsWrapperNode thread_node;
typedef serial_node NT;

typedef Tpetra::Operator<>::scalar_type SC;
typedef Tpetra::Operator<SC>::local_ordinal_type LO;
typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;
typedef Tpetra::Map<LO, GO, NT> map_type;
typedef Tpetra::Operator<SC, LO, GO, NT> op_type;
typedef Tpetra::MultiVector<SC, LO, GO, NT> mv_type;
typedef Tpetra::Vector<SC, LO, GO, NT> vec_type;
typedef Tpetra::CrsMatrix<SC, LO, GO, NT> matrix_type;
typedef Tpetra::CrsGraph<LO, GO, NT> graph_type;

typedef double ValueType;
typedef Kokkos::DefaultExecutionSpace DeviceSpaceType;

typedef Kokkos::DynRankView<ValueType, DeviceSpaceType> DynRankView;

typedef shards::CellTopology CellTopology;
// typedef Intrepid2::CellTools<DeviceSpaceType>    CellTools;
typedef Intrepid2::FunctionSpaceTools<DeviceSpaceType> fst;

using namespace std;
using namespace Intrepid2;

// Functions to evaluate exact solution and derivatives
KOKKOS_FUNCTION double evalu(double &x, double &y, double &z);
KOKKOS_FUNCTION int evalGradu(double &x, double &y, double &z, double &gradu1, double &gradu2, double &gradu3);
KOKKOS_FUNCTION double evalDivGradu(double &x, double &y, double &z);

int main(int argc, char *argv[]) {
    Tpetra::initialize(&argc, &argv);
    {
        // Check number of arguments
        if (argc < 4) {
            std::cout << "\n>>> ERROR: Invalid number of arguments.\n\n";
            std::cout << "Usage:\n\n";
            std::cout << "  ./poisson.exe NX NY NZ verbose\n\n";
            std::cout << " where \n";
            std::cout << "   int NX              - num intervals in x direction (assumed box domain, 0,1) \n";
            std::cout << "   int NY              - num intervals in y direction (assumed box domain, 0,1) \n";
            std::cout << "   int NZ              - num intervals in z direction (assumed box domain, 0,1) \n";
            std::cout << "   verbose (optional)  - any character, indicates verbose output \n\n";
            exit(1);
        }

        // This little trick lets us print to std::cout only if
        // a (dummy) command-line argument is provided.
        int iprint = argc - 1;
        Teuchos::RCP<std::ostream> outStream;
        Teuchos::oblackholestream bhs;  // outputs nothing
        if (iprint > 3)
            outStream = Teuchos::rcp(&std::cout, false);
        else
            outStream = Teuchos::rcp(&bhs, false);

        // Save the format state of the original std::cout.
        Teuchos::oblackholestream oldFormatState;
        oldFormatState.copyfmt(std::cout);

        *outStream << "===============================================================================\n"
                   << "|                                                                             |\n"
                   << "|  Example: Generate Stiffness Matrix and Right Hand Side Vector for          |\n"
                   << "|                   Poisson Equation on Hexahedral Mesh                       |\n"
                   << "|                                                                             |\n"
                   << "===============================================================================\n";

        // ************************************ GET INPUTS **************************************

        int NX = atoi(argv[1]);  // num intervals in x direction (assumed box domain, 0,1)
        int NY = atoi(argv[2]);  // num intervals in y direction (assumed box domain, 0,1)
        int NZ = atoi(argv[3]);  // num intervals in z direction (assumed box domain, 0,1)

        // *********************************** CELL TOPOLOGY **********************************

        // Get cell topology for base hexahedron
        CellTopology hex_8(shards::getCellTopologyData<shards::Hexahedron<8> >());

        // Get dimensions
        int numNodesPerElem = hex_8.getNodeCount();
        int spaceDim = hex_8.getDimension();

        // *********************************** GENERATE MESH ************************************

        *outStream << "Generating mesh ... \n\n";

        *outStream << "   NX"
                   << "   NY"
                   << "   NZ\n";
        *outStream << std::setw(5) << NX << std::setw(5) << NY << std::setw(5) << NZ << "\n\n";

        // Print mesh information
        int numElems = NX * NY * NZ;
        int numNodes = (NX + 1) * (NY + 1) * (NZ + 1);
        *outStream << " Number of Elements: " << numElems << " \n";
        *outStream << "    Number of Nodes: " << numNodes << " \n\n";

        // Cube
        double leftX = 0.0, rightX = 1.0;
        double leftY = 0.0, rightY = 1.0;
        double leftZ = 0.0, rightZ = 1.0;

        // Mesh spacing
        double hx = (rightX - leftX) / ((double)NX);
        double hy = (rightY - leftY) / ((double)NY);
        double hz = (rightZ - leftZ) / ((double)NZ);

        // Get nodal coordinates
        // FieldContainer<double> nodeCoord(numNodes, spaceDim);
        // FieldContainer<int> nodeOnBoundary(numNodes);
        DynRankView nodeCoord("nodeCoord", numNodes, spaceDim);
        Kokkos::DynRankView<int, DeviceSpaceType> nodeOnBoundary("nodeOnBoundary", numNodes);
        int inode = 0;
        for (int k = 0; k < NZ + 1; k++) {
            for (int j = 0; j < NY + 1; j++) {
                for (int i = 0; i < NX + 1; i++) {
                    nodeCoord(inode, 0) = leftX + (double)i * hx;
                    nodeCoord(inode, 1) = leftY + (double)j * hy;
                    nodeCoord(inode, 2) = leftZ + (double)k * hz;
                    if (k == 0 || j == 0 || i == 0 || k == NZ || j == NY || i == NX) {
                        nodeOnBoundary(inode) = 1;
                    } else {
                        nodeOnBoundary(inode) = 0;
                    }
                    inode++;
                }
            }
        }
#define DUMP_DATA
#ifdef DUMP_DATA
        // Print nodal coords
        ofstream fcoordout("coords.dat");
        for (int i = 0; i < numNodes; i++) {
            fcoordout << nodeCoord(i, 0) << " ";
            fcoordout << nodeCoord(i, 1) << " ";
            fcoordout << nodeCoord(i, 2) << "\n";
        }
        fcoordout.close();
#endif

        // Element to Node map
        // FieldContainer<int> elemToNode(numElems, numNodesPerElem);
        Kokkos::DynRankView<int, DeviceSpaceType> elemToNode("elemToNode", numElems, numNodesPerElem);
        int ielem = 0;
        for (int k = 0; k < NZ; k++) {
            for (int j = 0; j < NY; j++) {
                for (int i = 0; i < NX; i++) {
                    elemToNode(ielem, 0) = (NY + 1) * (NX + 1) * k + (NX + 1) * j + i;
                    elemToNode(ielem, 1) = (NY + 1) * (NX + 1) * k + (NX + 1) * j + i + 1;
                    elemToNode(ielem, 2) = (NY + 1) * (NX + 1) * k + (NX + 1) * (j + 1) + i + 1;
                    elemToNode(ielem, 3) = (NY + 1) * (NX + 1) * k + (NX + 1) * (j + 1) + i;
                    elemToNode(ielem, 4) = (NY + 1) * (NX + 1) * (k + 1) + (NX + 1) * j + i;
                    elemToNode(ielem, 5) = (NY + 1) * (NX + 1) * (k + 1) + (NX + 1) * j + i + 1;
                    elemToNode(ielem, 6) = (NY + 1) * (NX + 1) * (k + 1) + (NX + 1) * (j + 1) + i + 1;
                    elemToNode(ielem, 7) = (NY + 1) * (NX + 1) * (k + 1) + (NX + 1) * (j + 1) + i;
                    ielem++;
                }
            }
        }
#ifdef DUMP_DATA
        // Output connectivity
        ofstream fe2nout("elem2node.dat");
        for (int k = 0; k < NZ; k++) {
            for (int j = 0; j < NY; j++) {
                for (int i = 0; i < NX; i++) {
                    int ielem = i + j * NX + k * NX * NY;
                    for (int m = 0; m < numNodesPerElem; m++) {
                        fe2nout << elemToNode(ielem, m) << "  ";
                    }
                    fe2nout << "\n";
                }
            }
        }
        fe2nout.close();
#endif

        // ************************************ CUBATURE **************************************

        *outStream << "Getting cubature ... \n\n";

        // Get numerical integration points and weights
        DefaultCubatureFactory cubFactory;

        // auto cubature = DefaultCubatureFactory::create<DeviceSpaceType,ValueType,ValueType>(cellTopo, order);

        int cubDegree = 2;
        auto hexCub = cubFactory.create<DeviceSpaceType, ValueType, ValueType>(hex_8, cubDegree);

        int cubDim = hexCub->getDimension();
        int numCubPoints = hexCub->getNumPoints();

        // FieldContainer<double> cubPoints(numCubPoints, cubDim);
        // FieldContainer<double> cubWeights(numCubPoints);
        DynRankView cubPoints("cubPoints", numCubPoints, cubDim);
        DynRankView cubWeights("cubWeights", numCubPoints);

        hexCub->getCubature(cubPoints, cubWeights);

        // ************************************** BASIS ***************************************

        *outStream << "Getting basis ... \n\n";

        // Define basis
        Basis_HGRAD_HEX_C1_FEM<DeviceSpaceType, double, double> hexHGradBasis;

        //    Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > hexHGradBasis;
        int numFieldsG = hexHGradBasis.getCardinality();
        DynRankView hexGVals("hexGVals", numFieldsG, numCubPoints);
        DynRankView hexGrads("hexGrads", numFieldsG, numCubPoints, spaceDim);

        // Evaluate basis values and gradients at cubature points
        hexHGradBasis.getValues(hexGVals, cubPoints, OPERATOR_VALUE);
        hexHGradBasis.getValues(hexGrads, cubPoints, OPERATOR_GRAD);

        // ******** LOOP OVER ELEMENTS TO CREATE LOCAL STIFFNESS MATRIX *************

        *outStream << "Building stiffness matrix and right hand side ... \n\n";

        // Settings and data structures for mass and stiffness matrices

        // Container for nodes
        DynRankView hexNodes("hexNodes", numElems, numNodesPerElem, spaceDim);
        // Containers for Jacobian
        DynRankView hexJacobian("hexJacobian", numElems, numCubPoints, spaceDim, spaceDim);
        DynRankView hexJacobInv("hexJacobInv", numElems, numCubPoints, spaceDim, spaceDim);
        DynRankView hexJacobDet("hexJacobDet", numElems, numCubPoints);
        // Containers for element HGRAD stiffness matrix
        DynRankView localStiffMatrix("localStiffMatrix", numElems, numFieldsG, numFieldsG);
        DynRankView weightedMeasure("weightedMeasure", numElems, numCubPoints);
        DynRankView hexGradsTransformed("hexGradsTransformed", numElems, numFieldsG, numCubPoints, spaceDim);
        DynRankView hexGradsTransformedWeighted(
            "hexGradsTransformedWeighted", numElems, numFieldsG, numCubPoints, spaceDim);
        // Containers for right hand side vectors
        DynRankView rhsData("rhsData", numElems, numCubPoints);
        DynRankView localRHS("localRHS", numElems, numFieldsG);
        DynRankView hexGValsTransformed("hexGValsTransformed", numElems, numFieldsG, numCubPoints);
        DynRankView hexGValsTransformedWeighted("hexGValsTransformedWeighted", numElems, numFieldsG, numCubPoints);
        // Container for cubature points in physical space
        DynRankView physCubPoints("physCubPoints", numElems, numCubPoints, cubDim);
        int indexBase = 0;
        *outStream << "Parallel region ... \n\n";
        // Teuchos::GlobalMPISession mpiSession(&argc, &argv, nullptr);
        Teuchos::RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();
        Teuchos::RCP<const map_type> globalMapG = Teuchos::rcp(new map_type(numNodes, indexBase, Comm));
        Teuchos::RCP<matrix_type> StiffMatrix = Teuchos::rcp(new matrix_type(globalMapG, numFieldsG));
        Teuchos::RCP<vec_type> rhs = Teuchos::rcp(new vec_type(globalMapG, false));

        Kokkos::parallel_for(numElems, KOKKOS_LAMBDA(int k) {
            // *** Element loop ***
            // Physical cell coordinates
            for (int i = 0; i < numNodesPerElem; i++) {
                hexNodes(k, i, 0) = nodeCoord(elemToNode(k, i), 0);
                hexNodes(k, i, 1) = nodeCoord(elemToNode(k, i), 1);
                hexNodes(k, i, 2) = nodeCoord(elemToNode(k, i), 2);
            }
        });

        // Compute cell Jacobians, their inverses and their determinants
        Intrepid2::CellTools<DeviceSpaceType>::setJacobian(hexJacobian, cubPoints, hexNodes, hex_8);
        Intrepid2::CellTools<DeviceSpaceType>::setJacobianInv(hexJacobInv, hexJacobian);
        Intrepid2::CellTools<DeviceSpaceType>::setJacobianDet(hexJacobDet, hexJacobian);
        // ************************** Compute element HGrad stiffness matrices *******************************

        // transform to physical coordinates
        fst::HGRADtransformGRAD<double>(hexGradsTransformed, hexJacobInv, hexGrads);

        // compute weighted measure
        fst::computeCellMeasure<double>(weightedMeasure, hexJacobDet, cubWeights);

        // multiply values with weighted measure
        fst::multiplyMeasure<double>(hexGradsTransformedWeighted, weightedMeasure, hexGradsTransformed);

        // integrate to compute element stiffness matrix
        fst::integrate(localStiffMatrix, hexGradsTransformed, hexGradsTransformedWeighted);
        // assemble into global matrix

        // ******************************* Build right hand side ************************************

        // transform integration points to physical points
        Intrepid2::CellTools<DeviceSpaceType>::mapToPhysicalFrame(physCubPoints, cubPoints, hexNodes, hex_8);
        // evaluate right hand side function at physical points

        Kokkos::parallel_for(numElems, KOKKOS_LAMBDA(int k) {
            for (int nPt = 0; nPt < numCubPoints; nPt++) {
                double x = physCubPoints(k, nPt, 0);
                double y = physCubPoints(k, nPt, 1);
                double z = physCubPoints(k, nPt, 2);

                rhsData(k, nPt) = evalDivGradu(x, y, z);
            }
        });

        // transform basis values to physical coordinates
        fst::HGRADtransformVALUE<double>(hexGValsTransformed, hexGVals);

        // multiply values with weighted measure
        fst::multiplyMeasure<double>(hexGValsTransformedWeighted, weightedMeasure, hexGValsTransformed);
        // integrate rhs term
        fst::integrate(localRHS, rhsData, hexGValsTransformedWeighted);
        // assemble into global vector

        // auto local_mat = Teuchos::rcp(new KokkosSparse::CrsMatrix<double, int, DeviceSpaceType>());
        auto local_mat = StiffMatrix->getLocalMatrix();

        Kokkos::parallel_for(numElems, KOKKOS_LAMBDA(const int k) {
            ///
            for (int row = 0; row < numFieldsG; row++) {
                for (int col = 0; col < numFieldsG; col++) {
                    int rowIndex = elemToNode(k, row);
                    int colIndex = elemToNode(k, col);

                    double val = localStiffMatrix(k, row, col);

                    local_mat.sumIntoValues(rowIndex, &colIndex, 1, &val, true, true);
                }
            }
        });

        auto local_vec = rhs->getLocalView<DeviceSpaceType>();

        Kokkos::parallel_for(numElems, KOKKOS_LAMBDA(const int k) {
            for (int row = 0; row < numFieldsG; row++) {
                int rowIndex = elemToNode(k, row);
                const double val = -localRHS(k, row);

                // ATTENTION local index (second index because multivec)
                local_vec(rowIndex, 0) = val;
            }
        });

        // Assemble global matrices
        StiffMatrix->globalAssemble();
        StiffMatrix->fillComplete();
        //   rhs->globalAssemble();
        // Adjust stiffness matrix and rhs based on boundary conditions
        for (int row = 0; row < numNodes; row++) {
            if (nodeOnBoundary(row)) {
                int rowindex = row;
                for (int col = 0; col < numNodes; col++) {
                    double val = 0.0;
                    int colindex = col;
                    StiffMatrix->replaceGlobalValues(rowindex, 1, &val, &colindex);
                }
                double val = 1.0;

                StiffMatrix->replaceGlobalValues(rowindex, 1, &val, &rowindex);
                val = 0.0;
                rhs->replaceGlobalValue(rowindex, val);
            }
        }

        /*#ifdef DUMP_DATA
         // Dump matrices to disk
         EpetraExt::RowMatrixToMatlabFile("stiff_matrix.dat",StiffMatrix);
         EpetraExt::MultiVectorToMatrixMarketFile("rhs_vector.dat",rhs,0,0,false);
         #endif*/

        std::cout << "End Result: TEST PASSED\n";

        // reset format state of std::cout
        std::cout.copyfmt(oldFormatState);
    }
    Tpetra::finalize();
    return 0;
}

// Calculates value of exact solution u
KOKKOS_FUNCTION double evalu(double &x, double &y, double &z) {
    /*
     // function1
     double exactu = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
     */

    // function2
    double exactu = sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) * exp(x + y + z);

    return exactu;
}

// Calculates gradient of exact solution u
KOKKOS_FUNCTION int evalGradu(double &x, double &y, double &z, double &gradu1, double &gradu2, double &gradu3) {
    /*
     // function 1
     gradu1 = M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
     gradu2 = M_PI*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z);
     gradu3 = M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);
     */

    // function2
    gradu1 = (M_PI * cos(M_PI * x) + sin(M_PI * x)) * sin(M_PI * y) * sin(M_PI * z) * exp(x + y + z);
    gradu2 = (M_PI * cos(M_PI * y) + sin(M_PI * y)) * sin(M_PI * x) * sin(M_PI * z) * exp(x + y + z);
    gradu3 = (M_PI * cos(M_PI * z) + sin(M_PI * z)) * sin(M_PI * x) * sin(M_PI * y) * exp(x + y + z);

    return 0;
}

// Calculates Laplacian of exact solution u
KOKKOS_FUNCTION double evalDivGradu(double &x, double &y, double &z) {
    /*
     // function 1
     double divGradu = -3.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
     */

    // function 2
    double divGradu = -3.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) * exp(x + y + z) +
                      2.0 * M_PI * cos(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) * exp(x + y + z) +
                      2.0 * M_PI * cos(M_PI * y) * sin(M_PI * x) * sin(M_PI * z) * exp(x + y + z) +
                      2.0 * M_PI * cos(M_PI * z) * sin(M_PI * x) * sin(M_PI * y) * exp(x + y + z) +
                      3.0 * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) * exp(x + y + z);

    return divGradu;
}
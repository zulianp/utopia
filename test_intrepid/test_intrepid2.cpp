// @HEADER

/** \file   example_03.cpp
    \brief  Example building stiffness matrix and right hand side for a Poisson equation 
            using nodal (Hgrad) elements.

    \verbatim
             div grad u = f in Omega
                      u = 0 on Gamma 

     Discrete linear system for nodal coefficients(x):
        
                 Kx = b

            K - HGrad stiffness matrix
            b - right hand side vector 
                
    \endverbatim

    \author Created by P. Bochev, D. Ridzal and K. Peterson.

    
     \remark Usage
     \verbatim

     ./Intrepid_example_Drivers_Example_03.exe NX NY NZ verbose

        int NX              - num intervals in x direction (assumed box domain, 0,1)
        int NY              - num intervals in y direction (assumed box domain, 0,1)
        int NZ              - num intervals in z direction (assumed box domain, 0,1)
        verbose (optional)  - any character, indicates verbose output

     \endverbatim

    \remark Sample command line
    \code   ./Intrepid_example_Drivers_Example_03.exe 10 10 10 \endcode
*/

// Intrepid includes
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"

#include "Kokkos_Core.hpp"

// Tpetra includes
//#include "Tpetra_FECrsMatrix_decl.hpp"
//#include "Tpetra_FEVector_decl.hpp"
//#include "Tpetra_SerialComm.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_CrsMatrix.hpp>

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// TpetraExt includes
/*#include "TpetraExt_RowMatrixOut.hpp"
#include "TpetraExt_MultiVectorOut.hpp"*/

#include<chrono>

//using namespace std;
using namespace Intrepid2;

// Functions to evaluate exact solution and derivatives
double evalu(double & x, double & y, double & z);
int evalGradu(double & x, double & y, double & z, double & gradu1, double & gradu2, double & gradu3);
double evalDivGradu(double & x, double & y, double & z);

int main(int argc, char *argv[]) {


using Clock = std::chrono::high_resolution_clock; 
Clock::time_point start_time = Clock::now();

        //types of Operators
        typedef Tpetra::Operator<>::scalar_type SC;
        typedef Tpetra::Operator<SC>::local_ordinal_type LO;
        typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;


    //types of Kokkos Parallel Nodes
    typedef Kokkos::Compat::KokkosSerialWrapperNode serial_node;
    typedef Kokkos::Serial ks; //TODO difference between kokkos serial and serial node
#ifdef KOKKOS_ENABLE_CUDA
    typedef Kokkos::CUDA kc;
    typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
    typedef cuda_node NT;
#elif defined KOKKOS_ENABLE_ROCM //Kokkos::Compat::KokkosROCmWrapperNode doesn't exist
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::ROCm> rocm_node;
    typedef rocm_node NT;
#elif defined KOKKOS_ENABLE_OPENMP
    typedef Kokkos::OpenMP ko;
    typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
    typedef openmp_node NT;
#else
    typedef serial_node NT;
#endif


    typedef ks execution_space;

 Kokkos::initialize();
  //Check number of arguments
   if (argc < 4) {
      std::cout <<"\n>>> ERROR: Invalid number of arguments.\n\n";
      std::cout <<"Usage:\n\n";
      std::cout <<"  ./Intrepid_example_Drivers_Example_03.exe NX NY NZ verbose\n\n";
      std::cout <<" where \n";
      std::cout <<"   int NX              - num intervals in x direction (assumed box domain, 0,1) \n";
      std::cout <<"   int NY              - num intervals in y direction (assumed box domain, 0,1) \n";
      std::cout <<"   int NZ              - num intervals in z direction (assumed box domain, 0,1) \n";
      std::cout <<"   verbose (optional)  - any character, indicates verbose output \n\n";
      exit(1);
   }
  
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 3)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|  Example: Generate Stiffness Matrix and Right Hand Side Vector for          |\n" \
    << "|                   Poisson Equation on Hexahedral Mesh                       |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";


// ************************************ GET INPUTS **************************************

    int NX            = atoi(argv[1]);  // num intervals in x direction (assumed box domain, 0,1)
    int NY            = atoi(argv[2]);  // num intervals in y direction (assumed box domain, 0,1)
    int NZ            = atoi(argv[3]);  // num intervals in z direction (assumed box domain, 0,1)

// *********************************** CELL TOPOLOGY **********************************

   // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology hex_8(shards::getCellTopologyData<shards::Hexahedron<8> >() );

   // Get dimensions 
    int numNodesPerElem = hex_8.getNodeCount();
    int spaceDim = hex_8.getDimension();

// *********************************** GENERATE MESH ************************************

    *outStream << "Generating mesh ... \n\n";

    *outStream << "   NX" << "   NY" << "   NZ\n";
    *outStream << std::setw(5) << NX <<
                 std::setw(5) << NY <<
                 std::setw(5) << NZ << "\n\n";

   // Print mesh information
    int numElems = NX*NY*NZ;
    int numNodes = (NX+1)*(NY+1)*(NZ+1);
    *outStream << " Number of Elements: " << numElems << " \n";
    *outStream << "    Number of Nodes: " << numNodes << " \n\n";

   // Cube
    double leftX = 0.0, rightX = 1.0;
    double leftY = 0.0, rightY = 1.0;
    double leftZ = 0.0, rightZ = 1.0;

   // Mesh spacing
    double hx = (rightX-leftX)/((double)NX);
    double hy = (rightY-leftY)/((double)NY);
    double hz = (rightZ-leftZ)/((double)NZ);

   // Get nodal coordinates
    Kokkos::DynRankView<double,execution_space> nodeCoord("pota",numNodes, spaceDim);
    Kokkos::DynRankView<int,execution_space> nodeOnBoundary("pota1",numNodes);
    int inode = 0;
    for (int k=0; k<NZ+1; k++) {
      for (int j=0; j<NY+1; j++) {
        for (int i=0; i<NX+1; i++) {
          nodeCoord(inode,0) = leftX + (double)i*hx;
          nodeCoord(inode,1) = leftY + (double)j*hy;
          nodeCoord(inode,2) = leftZ + (double)k*hz;
          if (k==0 || j==0 || i==0 || k==NZ || j==NY || i==NX){
             nodeOnBoundary(inode)=1;
          }
          else {
             nodeOnBoundary(inode)=0;
          }
          inode++;
        }
      }
    }
#define DUMP_DATA
#ifdef DUMP_DATA
   // Print nodal coords
    std::ofstream fcoordout("coords.dat");
    for (int i=0; i<numNodes; i++) {
       fcoordout << nodeCoord(i,0) <<" ";
       fcoordout << nodeCoord(i,1) <<" ";
       fcoordout << nodeCoord(i,2) <<"\n";
    }
    fcoordout.close();
#endif


  // Element to Node map
    Kokkos::DynRankView<int,execution_space> elemToNode("elemToNode",numElems, numNodesPerElem);
    int ielem = 0;
    for (int k=0; k<NZ; k++) {
      for (int j=0; j<NY; j++) {
        for (int i=0; i<NX; i++) {
          elemToNode(ielem,0) = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i;
          elemToNode(ielem,1) = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i + 1;
          elemToNode(ielem,2) = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i + 1;
          elemToNode(ielem,3) = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i;
          elemToNode(ielem,4) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i;
          elemToNode(ielem,5) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i + 1;
          elemToNode(ielem,6) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i + 1;
          elemToNode(ielem,7) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i;
          ielem++;
        }
      }
    }
#ifdef DUMP_DATA
   // Output connectivity
    std::ofstream fe2nout("elem2node.dat");
    for (int k=0; k<NZ; k++) {
      for (int j=0; j<NY; j++) {
        for (int i=0; i<NX; i++) {
          int ielem = i + j * NX + k * NX * NY;
          for (int m=0; m<numNodesPerElem; m++){
              fe2nout << elemToNode(ielem,m) <<"  ";
           }
          fe2nout <<"\n";
        }
      }
    }
    fe2nout.close();
#endif


// ************************************ CUBATURE **************************************

    *outStream << "Getting cubature ... \n\n";

   // Get numerical integration points and weights
    DefaultCubatureFactory  cubFactory;                                   
    int cubDegree = 2;
    auto hexCub = cubFactory.create<execution_space,double,double>(hex_8, cubDegree); //Teuchos::RCP<Cubature<double> >

    int cubDim       = hexCub->getDimension();
    int numCubPoints = hexCub->getNumPoints();

    Kokkos::DynRankView<double,execution_space> cubPoints("cubPoints",numCubPoints, cubDim);
    Kokkos::DynRankView<double,execution_space> cubWeights("cubWeights",numCubPoints);
    hexCub->getCubature(cubPoints, cubWeights);


// ************************************** BASIS ***************************************

     *outStream << "Getting basis ... \n\n";

   // Define basis 
     Basis_HGRAD_HEX_C1_FEM<execution_space,double,double> hexHGradBasis;
     int numFieldsG = hexHGradBasis.getCardinality();
     Kokkos::DynRankView<double,execution_space> hexGVals("hexGVals",numFieldsG, numCubPoints); 
     Kokkos::DynRankView<double,execution_space> hexGrads("hexGrads",numFieldsG, numCubPoints, spaceDim); 

  // Evaluate basis values and gradients at cubature points
     hexHGradBasis.getValues(hexGVals, cubPoints, OPERATOR_VALUE);
     hexHGradBasis.getValues(hexGrads, cubPoints, OPERATOR_GRAD);


// ******** LOOP OVER ELEMENTS TO CREATE LOCAL STIFFNESS MATRIX *************

    *outStream << "Building stiffness matrix and right hand side ... \n\n";

 // Settings and data structures for mass and stiffness matrices
    typedef CellTools<execution_space>  CellTools;
    typedef FunctionSpaceTools<execution_space> fst;
    int numCells = 1; 

   // Container for nodes
    Kokkos::DynRankView<double,execution_space> hexNodes("hexNodes",numCells, numNodesPerElem, spaceDim);
   // Containers for Jacobian
    Kokkos::DynRankView<double,execution_space> hexJacobian("hexNodes",numCells, numCubPoints, spaceDim, spaceDim);
    Kokkos::DynRankView<double,execution_space> hexJacobInv("hexJacobInv",numCells, numCubPoints, spaceDim, spaceDim);
    Kokkos::DynRankView<double,execution_space> hexJacobDet("hexJacobDet",numCells, numCubPoints);
   // Containers for element HGRAD stiffness matrix
    Kokkos::DynRankView<double,execution_space> localStiffMatrix("localStiffMatrix",numCells, numFieldsG, numFieldsG);
    Kokkos::DynRankView<double,execution_space> weightedMeasure("weightedMeasure",numCells, numCubPoints);
    Kokkos::DynRankView<double,execution_space> hexGradsTransformed("hexGradsTransformed",numCells, numFieldsG, numCubPoints, spaceDim);
    Kokkos::DynRankView<double,execution_space> hexGradsTransformedWeighted("hexGradsTransformedWeighted",numCells, numFieldsG, numCubPoints, spaceDim);
   // Containers for right hand side vectors
    Kokkos::DynRankView<double,execution_space> rhsData("rhsData",numCells, numCubPoints);
    Kokkos::DynRankView<double,execution_space> localRHS("localRHS",numCells, numFieldsG);
    Kokkos::DynRankView<double,execution_space> hexGValsTransformed("hexGValsTransformed",numCells, numFieldsG, numCubPoints);
    Kokkos::DynRankView<double,execution_space> hexGValsTransformedWeighted("hexGValsTransformedWeighted",numCells, numFieldsG, numCubPoints);
   // Container for cubature points in physical space
    Kokkos::DynRankView<double,execution_space> physCubPoints("physCubPoints",numCells, numCubPoints, cubDim);

   // Global arrays in Tpetra format 
    Teuchos::RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm ();
    Teuchos::RCP<Tpetra::Map<LO, GO, NT> > globalMapG = Teuchos::rcp( new Tpetra::Map<LO, GO, NT>( numNodes, 0, Comm));

    Tpetra::CrsMatrix<SC, LO, GO, NT> StiffMatrix(globalMapG, numFieldsG);
    Tpetra::Vector<SC, LO, GO, NT> rhs(globalMapG, Teuchos::Copy);

 // *** Element loop ***
    for (int k=0; k<numElems; k++) {

     // Physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }

    // Compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobian, cubPoints, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInv, hexJacobian );
       CellTools::setJacobianDet(hexJacobDet, hexJacobian );

// ************************** Compute element HGrad stiffness matrices *******************************
  
     // transform to physical coordinates 
      fst::HGRADtransformGRAD<double>(hexGradsTransformed, hexJacobInv, hexGrads);
      
     // compute weighted measure
      fst::computeCellMeasure<double>(weightedMeasure, hexJacobDet, cubWeights);

     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGradsTransformedWeighted,
                                   weightedMeasure, hexGradsTransformed);

     // integrate to compute element stiffness matrix
      fst::integrate<double>(localStiffMatrix, hexGradsTransformed, hexGradsTransformedWeighted );

      // assemble into global matrix
      for (int row = 0; row < numFieldsG; row++){
        for (int col = 0; col < numFieldsG; col++){
            int rowIndex = elemToNode(k,row);
            int colIndex = elemToNode(k,col);
            double val = localStiffMatrix(0,row,col);
            StiffMatrix.insertGlobalValues(rowIndex, 1, &val, &colIndex);
         }
      }

// ******************************* Build right hand side ************************************

      // transform integration points to physical points
       CellTools::mapToPhysicalFrame(physCubPoints, cubPoints, hexNodes, hex_8);

      // evaluate right hand side function at physical points
       for (int nPt = 0; nPt < numCubPoints; nPt++){

          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);

          rhsData(0,nPt) = evalDivGradu(x, y, z);
       }

     // transform basis values to physical coordinates 
      fst::HGRADtransformVALUE<double>(hexGValsTransformed, hexGVals);

     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGValsTransformedWeighted,
                                   weightedMeasure, hexGValsTransformed);

     // integrate rhs term
      fst::integrate<double>(localRHS, rhsData, hexGValsTransformedWeighted );

    // assemble into global vector
     for (int row = 0; row < numFieldsG; row++){
           int rowIndex = elemToNode(k,row);
           double val = -localRHS(0,row);
           rhs.sumIntoGlobalValue(rowIndex, val);
      }
     
 } // *** end element loop ***


  // Assemble global matrices
   StiffMatrix.globalAssemble(); StiffMatrix.fillComplete();
  // rhs.globalAssemble();  //dove cazz.. è??  //TODO

 
  // Adjust stiffness matrix and rhs based on boundary conditions
   for (int row = 0; row<numNodes; row++){
       if (nodeOnBoundary(row)) {
          int rowindex = row;
          for (int col=0; col<numNodes; col++){
              double val = 0.0;
              int colindex = col;
              StiffMatrix.replaceGlobalValues(rowindex, 1, &val, &colindex);
          }
          double val = 1.0;
          StiffMatrix.replaceGlobalValues(rowindex, 1, &val, &rowindex);
          val = 0.0;
          rhs.replaceGlobalValue(rowindex, val);
       }
    }

#ifdef DUMP_DATA
   // Dump matrices to disk //TODO
 //    TpetraExt::RowMatrixToMatlabFile("stiff_matrix.dat",StiffMatrix);
 //    TpetraExt::MultiVectorToMatrixMarketFile("rhs_vector.dat",rhs,0,0,false);
#endif

   std::cout << "End Result: TEST PASSED\n";
   double time_spent = std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - start_time).count();
   std::cout << "Total Time: " << time_spent << std::endl;
   
   // reset format state of std::cout
   std::cout.copyfmt(oldFormatState);
   return 0;
   
   Kokkos::finalize();
}


// Calculates value of exact solution u
 double evalu(double & x, double & y, double & z)
 {
 /*
   // function1
    double exactu = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
 */

   // function2
   double exactu = std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z)*std::exp(x+y+z);

   return exactu;
 }

// Calculates gradient of exact solution u
 int evalGradu(double & x, double & y, double & z, double & gradu1, double & gradu2, double & gradu3)
 {
 /*
   // function 1
       gradu1 = M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
       gradu2 = M_PI*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z);
       gradu3 = M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);
 */

   // function2
       gradu1 = (M_PI*std::cos(M_PI*x)+std::sin(M_PI*x))
                  *std::sin(M_PI*y)*std::sin(M_PI*z)*std::exp(x+y+z);
       gradu2 = (M_PI*std::cos(M_PI*y)+std::sin(M_PI*y))
                  *std::sin(M_PI*x)*std::sin(M_PI*z)*std::exp(x+y+z);
       gradu3 = (M_PI*std::cos(M_PI*z)+std::sin(M_PI*z))
                  *std::sin(M_PI*x)*std::sin(M_PI*y)*std::exp(x+y+z);
  
   return 0;
 }

// Calculates Laplacian of exact solution u
 double evalDivGradu(double & x, double & y, double & z)
 {
 /*
   // function 1
    double divGradu = -3.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
 */

   // function 2
   double divGradu = -3.0*M_PI*M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z)*exp(x+y+z)
                    + 2.0*M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z)*std::exp(x+y+z)
                    + 2.0*M_PI*std::cos(M_PI*y)*std::sin(M_PI*x)*std::sin(M_PI*z)*std::exp(x+y+z)
                    + 2.0*M_PI*std::cos(M_PI*z)*std::sin(M_PI*x)*std::sin(M_PI*y)*std::exp(x+y+z)
                    + 3.0*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z)*std::exp(x+y+z);
   
   return divGradu;
 }


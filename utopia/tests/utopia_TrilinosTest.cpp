#include "utopia_Base.hpp"
//#include "utopia_trilinos_Vector.hpp"
#ifdef WITH_TRILINOS
#include "utopia_TrilinosTest.hpp"
#include "utopia.hpp"

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

  if (comm->getRank () == 0) {
    // On (MPI) Process 0, print out the Tpetra software version.
    std::cout << Tpetra::version () << std::endl << std::endl;
  }  
  
  // create the Map
  Teuchos::RCP<const map_type> map = rcp (new map_type (numGblIndices, indexBase, comm));
  const size_t numMyElements = map->getNodeNumElements ();

    
  std::cout<< "Creating the sparse matrix" << std::endl;
  // Create a Tpetra sparse matrix whose rows have distribution given by the Map.
  Teuchos::RCP<crs_matrix_type> A (new crs_matrix_type (map, 0));

    
  // Fill the sparse matrix, one row at a time.
  const scalar_type two    = static_cast<scalar_type> (2.0);
  const scalar_type negOne = static_cast<scalar_type> (-1.0);
  Teuchos::RCP<const Teuchos::Comm<int> > comm2 = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));


//  utopia::TpetraMatrix A_void2();
  utopia::TpetraMatrix A_void(comm2);
  utopia::TpetraMatrix A_map(map);
  utopia::TpetraMatrix A_mat(A_map);

  utopia::TpetraVector x_void();
  utopia::TpetraVector x_map(map);

//size_t numMyElements_void = A_void2.getNodeNumElements();

 for (local_ordinal_type lclRow = 0; lclRow < static_cast<local_ordinal_type> (numMyElements); ++lclRow) { 
    const global_ordinal_type gblRow = map->getGlobalElement (lclRow);
     A->insertGlobalValues (gblRow, Teuchos::tuple<global_ordinal_type> (gblRow, gblRow + 1), Teuchos::tuple<scalar_type> (two, negOne)); 
     A_map.insertGlobalValues (gblRow, Teuchos::tuple<global_ordinal_type> (gblRow, gblRow + 1), Teuchos::tuple<scalar_type> (two, negOne));
  }



const size_t numMyElements_void = A_void.getNodeNumElements();
 for (local_ordinal_type lclRow = 0; lclRow < static_cast<local_ordinal_type> (numMyElements_void); ++lclRow) { 
    const global_ordinal_type gblRow = A_void.getGlobalElement (lclRow);
    A_void.insertGlobalValues (gblRow, Teuchos::tuple<global_ordinal_type> (gblRow, gblRow + 1), Teuchos::tuple<scalar_type> (two, negOne));
  }

  // Tell the sparse matrix that we are done adding entries to it.
  A->fillComplete();
  A_void.fillComplete();
  A_map.fillComplete();


    // Run the power method and report the result.
  scalar_type lambda;// = PowerMethod<crs_matrix_type>::run (*A_map, niters, tolerance, out);

  std::cout << "Estimated max eigenvalue: " << lambda << std::endl;
  

}

namespace utopia {
	void run_trilinos_test()
	{
		UTOPIA_UNIT_TEST_BEGIN("TrilinosTest");
		UTOPIA_RUN_TEST(trilinos_example_test);
		UTOPIA_UNIT_TEST_END("TrilinosTest");
	}
}

#else //WITH_TRILINOS
namespace utopia {
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

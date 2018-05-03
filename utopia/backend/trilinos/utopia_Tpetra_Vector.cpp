#include "utopia_Tpetra_Vector.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Instance.hpp"

#include <Tpetra_CrsMatrix_decl.hpp>
#include <MatrixMarket_Tpetra.hpp>

namespace utopia {

	// void TpetraVector::values(const TpetraVector::rcp_comm_type &comm, std::size_t n_local, Tpetra::global_size_t n_global, TpetraVector::Scalar value)
	// {
	//     rcp_map_type map;

	//     if(n_local == INVALID_INDEX) {
	//         map = Teuchos::rcp(new map_type(n_global, 0, comm));
	//     } else {
	//         map = Teuchos::rcp(new map_type(n_global, n_local, 0, comm));
	//     }

	//     vec_ = Teuchos::rcp(new vector_type(map));
	//     implementation().putScalar(value);
	// }


	bool TpetraVector::read(const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const std::string &path)
	{
		typedef Tpetra::CrsMatrix<>                       crs_matrix_type;

		try {
			rcp_map_type map;
			vec_ = Tpetra::MatrixMarket::Reader<crs_matrix_type>::readVectorFile(path, comm, map);
		} catch(const std::exception &ex) {
			// is.close();
			std::cout << ex.what() << std::endl;
			return false;
		}

		return !vec_.is_null();
	}

	bool TpetraVector::write(const std::string &path) const
	{
		typedef Tpetra::CrsMatrix<>                       crs_matrix_type;
		
		if(vec_.is_null()) return false;

		try {
			Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeDenseFile(path, vec_, "vec", "");
		} catch(const std::exception &ex) {
			std::cout << ex.what() << std::endl;
			return false;
		}

		return true;
	}
}

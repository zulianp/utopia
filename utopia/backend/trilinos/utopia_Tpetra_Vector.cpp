#include "utopia_Tpetra_Vector.hpp"

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
		// Tpetra::MatrixMarket::Reader< decltype(m.implementation()) >
		return false;
	}

	bool TpetraVector::write(const std::string &path) const
	{
		// Tpetra::MatrixMarket::Reader< decltype(m.implementation()) >

		return false;
	}
}

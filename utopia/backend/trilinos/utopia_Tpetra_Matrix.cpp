#include "utopia_Tpetra_Matrix.hpp"

namespace utopia {
	void TpetraMatrix::mult(const TpetraVector &vec, TpetraVector &result) const
	{
		if(result.is_null()) {
			result.init(mat_->getRowMap());
		} else if(!result.implementation().getMap()->isSameAs(*mat_->getRowMap())) {
			result.init(mat_->getRowMap());
		}

		mat_->apply(vec.implementation(), result.implementation());
	}

	void TpetraMatrix::crs_init(
	              const rcp_comm_type &comm,
	              std::size_t rows_local,
	              std::size_t cols_local,
	              Tpetra::global_size_t rows_global,
	              Tpetra::global_size_t cols_global,
	              std::size_t nnz_x_row)
	{
		//Trilinos has more distribution options than petsc and it does not require to have 
		//a column operator structure as petsc
		if(cols_global == Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()) {
			Tpetra::global_size_t send_buff = cols_local;
			cols_global = 0;
			Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &send_buff, &cols_global);
		}

	    owner_ = true;
	    auto row_map = Teuchos::rcp(new map_type(rows_global, rows_local, 0, comm));
	    auto col_map = Teuchos::rcp(new map_type(cols_global, 0, comm, Tpetra::LocallyReplicated));
	    mat_ = Teuchos::rcp(new crs_matrix_type(row_map, nnz_x_row, Tpetra::DynamicProfile));
	    mat_->replaceColMap(col_map);
	}
}
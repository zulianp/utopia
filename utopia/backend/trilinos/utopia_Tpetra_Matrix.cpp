#include "utopia_Tpetra_Matrix.hpp"
#include <TpetraExt_MatrixMatrix_def.hpp>
#include "utopia_Logger.hpp"
#include "utopia_Instance.hpp"

namespace utopia {

	void TpetraMatrix::set(const global_ordinal_type &row, const global_ordinal_type &col, const Scalar &value)
	{
	    m_utopia_status_once(
	    	"> TpetraMatrix::set does what is supposed to do with respect to the edsl. " 
	    	"However it might be slow because of the double trilinos call and branching."
	    );

	    if(implementation().replaceGlobalValues(row, 1, &value, &col) != 1) {
	        implementation().insertGlobalValues(row, 1, &value, &col);
	    }
	}
	
	void TpetraMatrix::add(const global_ordinal_type &row, const global_ordinal_type &col, const Scalar &value)
	{
		m_utopia_status_once(
			"> TpetraMatrix::add does what is supposed to do with respect to the edsl. " 
			"However it might be slow because of the double trilinos call and branching."
		);

	    if(implementation().sumIntoGlobalValues(row, 1, &value, &col, false) != 1) {
	        implementation().insertGlobalValues(row, 1, &value, &col);
	    }
	}

	void TpetraMatrix::mult(const TpetraVector &vec, TpetraVector &result) const
	{
		if(result.is_null()) {
			result.init(mat_->getRowMap());
		} else if(!result.implementation().getMap()->isSameAs(*mat_->getRowMap())) {
			result.init(mat_->getRowMap());
		}

		mat_->apply(vec.implementation(), result.implementation());
	}

	void TpetraMatrix::mult(const TpetraMatrix &right, TpetraMatrix &result) const
	{
		if(result.is_null()) {
			auto s = size();
			auto ls = local_size();
			// auto row_map = Teuchos::rcp(new map_type(s.get(0), ls.get(0), 0, communicator()));
			// result.mat_  = Teuchos::rcp(new crs_matrix_type(row_map, 0, Tpetra::DynamicProfile));
			result.mat_ = Teuchos::rcp(new crs_matrix_type(implementation().getRowMap(), 0, Tpetra::DynamicProfile));
		}

		Tpetra::MatrixMatrix::Multiply(
			this->implementation(),
			false,
			right.implementation(),
			false,
			result.implementation()
		);
	}

	void TpetraMatrix::axpy(const Scalar alpha, const TpetraMatrix &x)
	{
		// write_lock();
		
		// Add (const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, bool transposeA, Scalar scalarA, CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, Scalar scalarB)
		
		// write_lock();
		//THIS DOES NOT WORK??????
		// Tpetra::MatrixMatrix::Add(
		// 	x.implementation(),
		// 	false,
		// 	alpha,
		// 	implementation(),
		// 	1.
		// );

		// write_unlock();

		Tpetra::MatrixMatrix::add(
			alpha,
			false,
			x.implementation(),
			1.,
			false,
			implementation()
		);

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

		rcp_map_type row_map;

		if(rows_local == INVALID_INDEX) {
			row_map = Teuchos::rcp(new map_type(rows_global, 0, comm));
		} else {
			if(cols_global == Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()) {
				Tpetra::global_size_t send_buff = cols_local;
				cols_global = 0;
				Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &send_buff, &cols_global);
			}

			row_map = Teuchos::rcp(new map_type(rows_global, rows_local, 0, comm));
		}

	    auto col_map = Teuchos::rcp(new map_type(cols_global, 0, comm, Tpetra::LocallyReplicated));
	    mat_ = Teuchos::rcp(new crs_matrix_type(row_map, nnz_x_row, Tpetra::DynamicProfile));
	    mat_->replaceColMap(col_map);

	    owner_ = true;
	}
}

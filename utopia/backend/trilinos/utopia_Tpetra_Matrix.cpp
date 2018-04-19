#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Instance.hpp"

#include <TpetraExt_MatrixMatrix_def.hpp>
#include <Tpetra_RowMatrixTransposer_decl.hpp>
#include <MatrixMarket_Tpetra.hpp>

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
			// result.owner_ = true;
		} else if(!result.implementation().getMap()->isSameAs(*mat_->getRowMap())) {
			result.init(mat_->getRowMap());
			// result.owner_ = true;
		}

		mat_->apply(vec.implementation(), result.implementation());
	}

	void TpetraMatrix::mult(const TpetraMatrix &right, TpetraMatrix &result) const
	{
		mult(false, right, false, result);
	}

	void TpetraMatrix::mult_t(const TpetraMatrix &right, TpetraMatrix &result) const
	{
		mult(true, right, false, result);
	}

	//result op(*this) * op
	void TpetraMatrix::mult(const bool transpose_this, const TpetraMatrix &right, const bool transpose_right, TpetraMatrix &result) const
	{
		if(result.is_null()) {
			auto s = size();
			auto ls = local_size();
			// auto row_map = Teuchos::rcp(new map_type(s.get(0), ls.get(0), 0, communicator()));
			// result.mat_  = Teuchos::rcp(new crs_matrix_type(row_map, 0, Tpetra::DynamicProfile));
			result.mat_ = Teuchos::rcp(new crs_matrix_type(implementation().getRowMap(), 0, Tpetra::DynamicProfile));
			result.owner_ = true;
		}

		//C = op(A)*op(B), 
		Tpetra::MatrixMatrix::Multiply(
			this->implementation(),
			transpose_this,
			right.implementation(),
			transpose_right,
			result.implementation()
		);
	}

	void TpetraMatrix::transpose(TpetraMatrix &mat) const
	{
		Tpetra::RowMatrixTransposer<Scalar, local_ordinal_type, global_ordinal_type, node_type> transposer(mat_);
		mat.mat_ = transposer.createTranspose();
		mat.owner_ = true;
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

	void TpetraMatrix::get_diag(TpetraVector &d) const
	{
		if(d.is_null()) {
			m_utopia_warning_once("TpetraMatrix::get_diag Assuming row <= col");
			d.init(implementation().getRowMap());
		}

		implementation().getLocalDiagCopy(d.implementation());
	}

	void TpetraMatrix::init_diag(const TpetraVector &d)
	{
		//FIXME maybe there is a better and more efficent way to do this
		//also without const_cast

		auto ls = d.local_size().get(0);
		auto gs = d.size().get(0);

		crs_init(d.communicator(),
				 ls,
				 ls,
				 gs,
				 gs,
				 1);


		auto r = d.range();

		const_cast<TpetraVector &>(d).read_lock();
		write_lock();

		for(auto i = r.begin(); i < r.end(); ++i) {
			set(i, i, d.get(i));
		}

		const_cast<TpetraVector &>(d).read_unlock();
		write_unlock();
	}

	bool TpetraMatrix::read(const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const std::string &path)
	{
		std::ifstream is;
		is.open(path.c_str());

		if(!is.good()) {
			return false;
		}
		
		try {
			//https://people.sc.fsu.edu/~jburkardt/data/mm/mm.html
			mat_ = Tpetra::MatrixMarket::Reader<crs_matrix_type>::readSparse(is, comm);
		} catch(std::exception &ex) {
			is.close();
			std::cout << ex.what() << std::endl;
			return false;
		}

		is.close();
		return !mat_.is_null();
	}

	bool TpetraMatrix::write(const std::string &path) const
	{
		// Tpetra::MatrixMarket::Reader< decltype(m.implementation()) >

		return false;
	}
}

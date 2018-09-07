#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Instance.hpp"

#include <TpetraExt_MatrixMatrix_def.hpp>
#include <Tpetra_RowMatrixTransposer_decl.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <iterator>

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


	//FIXME make faster version by storing view?
	TpetraMatrix::Scalar TpetraMatrix::get(const global_ordinal_type &row, const global_ordinal_type &col) const
	{
		Teuchos::ArrayView<const global_ordinal_type> cols;
		Teuchos::ArrayView<const Scalar> values;

		assert(implementation().isLocallyIndexed());

		auto local_col = col - implementation().getColMap()->getMinGlobalIndex();

		auto rr = row_range();
		implementation().getLocalRowView(row - rr.begin(), cols, values);

		auto it = std::lower_bound(std::begin(cols), std::end(cols), local_col);

		if(it == std::end(cols)) {
			return 0.;
		}

		assert(it != std::end(cols));

		std::size_t index = std::distance(std::begin(cols), it);

		assert(cols[index] == local_col);

		return values[index];
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
		try {
			mat_->apply(vec.implementation(), result.implementation());
		} catch(const std::exception &ex) {
			std::cout << ex.what() << std::endl;
			assert(false);
		}
	}

	void TpetraMatrix::mult_t(const TpetraVector &vec, TpetraVector &result) const
	{
		assert(mat_->hasTransposeApply());

		if(result.is_null()) {
			result.init(mat_->getDomainMap());
			// result.owner_ = true;
		} else if(!result.implementation().getMap()->isSameAs(*mat_->getDomainMap())) {
			result.init(mat_->getDomainMap());
			// result.owner_ = true;
		}
		try {
			mat_->apply(vec.implementation(), result.implementation(), Teuchos::TRANS);
		} catch(const std::exception &ex) {
			std::cout << ex.what() << std::endl;
			assert(false);
		}
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
		m_utopia_status_once("TpetraMatrix::mult Proper thing to do would be to check if the maps are compatible");
		//IMPROVEME
		result.mat_.reset();

		assert(!transpose_right);

		assert(transpose_this  || (this->local_size().get(1) == right.local_size().get(0) && this->size().get(1) == right.size().get(0) ) );
		assert(!transpose_this || (this->local_size().get(0) == right.local_size().get(0) && this->size().get(0) == right.size().get(0) ) );

		// if(result.is_null()) {

		auto col_map = Teuchos::rcp(new map_type(right.size().get(1), 0, communicator(), Tpetra::LocallyReplicated));

		if(transpose_this) {

			assert(!right.implementation().getDomainMap().is_null());
			result.mat_ = Teuchos::rcp(
				new crs_matrix_type(
					implementation().getDomainMap(),
					col_map,
					0, Tpetra::DynamicProfile));

		} else {

			result.mat_ = Teuchos::rcp(
				new crs_matrix_type(
					implementation().getRowMap(),
					col_map,
					0, Tpetra::DynamicProfile));
		}

		result.owner_ = true;
		// }

		try {
				//C = op(A)*op(B),
				Tpetra::MatrixMatrix::Multiply(
					this->implementation(),
					transpose_this,
					right.implementation(),
					transpose_right,
					result.implementation(),
					false
				);


			auto dm = this->implementation().getDomainMap();
			auto rm = this->implementation().getRangeMap();

			result.implementation().fillComplete(
				right.implementation().getDomainMap(),
				(transpose_this ? dm : rm)
			);

			// std::cout << ("---------------------------") << std::endl;
			// std::cout << transpose_this << std::endl;
			// disp(this->size());
			// disp(right.size());
			// disp(result.size());
			// std::cout << ("---------------------------") << std::endl;
			// disp(this->local_size());
			// disp(right.local_size());
			// disp(result.local_size());
			// std::cout << ("---------------------------") << std::endl;

			assert(transpose_this  || (this->local_size().get(0) == result.local_size().get(0) && this->size().get(0) == result.size().get(0) ) );
			assert(!transpose_this || (this->local_size().get(1) == result.local_size().get(0) && this->size().get(1) == result.size().get(0) ) );
			assert(transpose_right || (right.local_size().get(1) == result.local_size().get(1) && right.size().get(1) == result.size().get(1) ) );


		} catch(const std::exception &ex) {
			std::cout << ex.what() << std::endl;
			assert(false);
		}
	}

	void TpetraMatrix::transpose(TpetraMatrix &mat) const
	{
		//FIXME this does not work as it should
		try {
			Tpetra::RowMatrixTransposer<Scalar, local_ordinal_type, global_ordinal_type, node_type> transposer(mat_);
			mat.mat_ = transposer.createTranspose();
			mat.owner_ = true;


			//None of this creat a valid matrix for getGlobalRowView
			//1)
			// auto col_map = Teuchos::rcp(new map_type(size().get(0), 0, communicator(), Tpetra::LocallyReplicated));
			// mat.mat_->replaceColMap(col_map);

			//2)
			// mat.implementation().resumeFill();
			// mat.implementation().fillComplete(this->implementation().getRangeMap(), this->implementation().getDomainMap());

			assert(this->local_size().get(0) == mat.local_size().get(1));
			assert(this->local_size().get(1) == mat.local_size().get(0));

			assert(this->size().get(0) == mat.size().get(1));
			assert(this->size().get(1) == mat.size().get(0));

			assert(is_valid(true));
		} catch(const std::exception &ex) {
			std::cout << ex.what() << std::endl;
			assert(false);
		}
	}

	void TpetraMatrix::axpy(const Scalar alpha, const TpetraMatrix &x)
	{
		try {
			mat_ = Tpetra::MatrixMatrix::add(
				alpha,
				false,
				x.implementation(),
				1.,
				false,
				implementation()
			);

			owner_ = true;
		} catch(const std::exception &ex) {
			std::cout << ex.what() << std::endl;
			assert(false);
		}
	}

	void TpetraMatrix::finalize()
    {
    	try {
	    	if(init_) {

	    		assert(!init_->domain_map.is_null());
	    		assert(!init_->range_map.is_null());

	    		implementation().fillComplete(init_->domain_map, init_->range_map);
	    		// init_.reset();
	    	} else {
	    		// assert(false);
	        	implementation().fillComplete(implementation().getDomainMap(), implementation().getRangeMap());
	        }
        } catch(const std::exception &ex) {
        	std::cout << ex.what() << std::endl;
        	assert(false);
        }
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
			row_map = Teuchos::rcp(new map_type(rows_global, rows_local, 0, comm));
		}

		if(cols_global == Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()) {
			Tpetra::global_size_t send_buff = cols_local;
			cols_global = 0;
			Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &send_buff, &cols_global);
		}

	    auto col_map = Teuchos::rcp(new map_type(cols_global, 0, comm, Tpetra::LocallyReplicated));
	    mat_ = Teuchos::rcp(new crs_matrix_type(row_map, col_map, nnz_x_row, Tpetra::DynamicProfile));
	    owner_ = true;

	    init_ = std::make_shared<InitStructs>();
	    if(cols_local == INVALID_INDEX) {
	    	init_->domain_map = Teuchos::rcp(new map_type(cols_global, 0, comm));
	    } else {
	    	init_->domain_map = Teuchos::rcp(new map_type(cols_global, cols_local, 0, comm));
	    }

	    init_->range_map = row_map;
	}

	void TpetraMatrix::crs_identity(const rcp_comm_type &comm,
	              std::size_t rows_local,
	              std::size_t cols_local,
	              Tpetra::global_size_t rows_global,
	              Tpetra::global_size_t cols_global,
	              const Scalar factor)
	{
		crs_init(comm, rows_local, cols_local, rows_global, cols_global, 1.);

		write_lock();

		Range r = row_range();

		auto cols = init_->domain_map->getGlobalNumElements();

		for(auto i = r.begin(); i < r.end(); ++i) {
			if(i < cols) {
				set(i, i, factor);
			} else {
				break;
			}
		}

		write_unlock();
	}

	void TpetraMatrix::get_diag(TpetraVector &d) const
	{
		const bool is_row_min = this->size().get(0) <= this->size().get(1);
		global_ordinal_type n = (is_row_min)? this->size().get(0) : this->size().get(1);

		if(d.is_null() || d.size().get(0) != n) {
			m_utopia_warning_once("TpetraMatrix::get_diag Assuming row <= col");

			if(is_row_min) {
				d.init(implementation().getRowMap());
			} else {
				d.init(implementation().getDomainMap());
			}
		}

		implementation().getLocalDiagCopy(d.implementation());
	}

	void TpetraMatrix::init_diag(const TpetraVector &d)
	{
		auto ls = d.local_size().get(0);
		auto gs = d.size().get(0);

		crs_init(d.communicator(),
				 ls,
				 ls,
				 gs,
				 gs,
				 1);


		auto r = d.range();
		auto data = d.implementation().getData();

		assert(!data.is_null());

		write_lock();

		local_ordinal_type index = 0;

		for(auto i = r.begin(); i < r.end(); ++i) {
			set(i, i, data[index++]);
		}

		write_unlock();
	}

	bool TpetraMatrix::read(const Teuchos::RCP< const Teuchos::Comm<int> > &comm, const std::string &path)
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
		if(mat_.is_null()) return false;

		try {
			Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparseFile(path, mat_, "mat", "", false);
		} catch(const std::exception &ex) {
			std::cout << ex.what() << std::endl;
			return false;
		}

		return true;
	}

	bool TpetraMatrix::is_valid(const bool verbose) const
	{
		if(mat_.is_null()) {
			if(verbose) { std::cerr << "is_null" << std::endl; }
			return false;
		}

		auto comm = communicator();

		if(comm->getSize() == 1) {

			if(local_size() != size()) {

				if(verbose) {
					std::cerr << "local_size() != size()" << std::endl;
					std::cerr << local_size() << " != " << size() << std::endl;
					std::cerr << "this indicates inconsistent domain_map with respect to the col_map" << std::endl;
				}

				return false;
			}
		}

		return true;
	}

	TpetraMatrix::Scalar TpetraMatrix::norm2() const
	{
		return implementation().getFrobeniusNorm();
	}

}

#include "utopia_Tpetra_Vector.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Instance.hpp"
#include "utopia_kokkos_Eval_Reduce.hpp"

#include <Tpetra_CrsMatrix_decl.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <Kokkos_Core.hpp>

#include <cmath>

namespace utopia {

	void TpetraVector::add_vector(
	    const std::vector<GO> &indices,
	    const std::vector<Scalar> &values)
	{
		const std::size_t n = values.size();
		assert(n == indices.size());

		for(std::size_t i = 0; i < n; ++i) {
			add(indices[i], values[i]);
		}
	}

	void TpetraVector::set_vector(
	    const std::vector<GO> &indices,
	    const std::vector<Scalar> &values)
	{
		const std::size_t n = values.size();
		assert(n == indices.size());

		for(std::size_t i = 0; i < n; ++i) {
			set(indices[i], values[i]);
		}
	}

	bool TpetraVector::read(const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const std::string &path)
	{
        typedef Tpetra::Operator<>::scalar_type SC;
        typedef Tpetra::Operator<SC>::local_ordinal_type LO;
        typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;

        typedef Tpetra::CrsMatrix<SC, LO, GO, NT>         crs_matrix_type;

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
     typedef Tpetra::CrsMatrix<SC, LO, GO, NT>            crs_matrix_type;
		if(vec_.is_null()) return false;

		try {
			Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeDenseFile(path, vec_, "vec", "");
		} catch(const std::exception &ex) {
			std::cout << ex.what() << std::endl;
			return false;
		}

		return true;
	}

	TpetraVector::Scalar TpetraVector::sum() const
	{
	    Scalar ret = KokkosEvalReduce<TpetraVector, Plus>::eval(*this, Plus(), Scalar(0.));
	    auto &comm = *communicator();
	    Scalar ret_global = 0.;
	    Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 1, &ret, &ret_global);
	    return ret_global;
	}

	TpetraVector::Scalar TpetraVector::min() const
	{
	    Scalar ret = KokkosEvalReduce<TpetraVector, Min>::eval(*this, Min(), std::numeric_limits<Scalar>::max());
	    auto &comm = *communicator();
	    Scalar ret_global = 0.;
	    Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, 1, &ret, &ret_global);
	    return ret_global;
	}

	TpetraVector::Scalar TpetraVector::max() const
	{
	  	Scalar ret = KokkosEvalReduce<TpetraVector, Max>::eval(*this, Max(), -std::numeric_limits<Scalar>::max());
	  	auto &comm = *communicator();
	  	Scalar ret_global = 0.;
	    Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &ret, &ret_global);
	    return ret_global;
	}

	bool TpetraVector::is_nan_or_inf() const
	{
		int ret = KokkosEvalReduce<TpetraVector, IsNaNOrInf>::eval(*this, IsNaNOrInf(), Scalar(0));
		auto &comm = *communicator();
		int ret_global = 0;

		Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &ret, &ret_global);
		return ret_global;
	}

	void TpetraVector::ghosted(
		const rcp_comm_type &comm, 
		const TpetraVector::GO &local_size,
	    const TpetraVector::GO &global_size,
		const std::vector<GO> &ghost_index
	)
	{
		rcp_map_type map(new map_type(global_size, local_size, 0, comm));
		rcp_map_type ghost_map;

		if(!ghost_index.empty()) {

			Range r = { map->getMinGlobalIndex(), map->getMaxGlobalIndex() + 1 };

			std::vector<GO> filled_with_local;
			filled_with_local.reserve(r.extent() + ghost_index.size());

			for(auto i = r.begin(); i != r.end(); ++i) {
				filled_with_local.push_back(i);
			}

			for(auto g : ghost_index) {
				if(!r.inside(g)) {
					filled_with_local.push_back(g);
				}
			}

			const Teuchos::ArrayView<const GO>
			   local_indices(filled_with_local);

			 GO local_entries = local_indices.size();
			 GO total_entries = 0;

			 Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &local_entries, &total_entries);

			 ghost_map = Teuchos::rcp(new map_type(total_entries, local_indices, 0, comm));

		} else {
			ghost_map = map;
		}

		ghosted_vec_ = Teuchos::rcp(new vector_type(ghost_map, 1));
		vec_ = ghosted_vec_->offsetViewNonConst(map, 0);

		assert(!vec_.is_null());
	}

	void TpetraVector::update_ghosts()
	{
		if(!has_ghosts()) return;

		auto map = vec_->getMap();
		auto ghost_map = ghosted_vec_->getMap();

		Tpetra::Import<
			LO,
			GO,
			vector_type::node_type> importer(map, ghost_map);

		ghosted_vec_->doImport(*vec_, importer, Tpetra::INSERT);
	}

	void TpetraVector::export_ghosts_add()
	{
		if(!has_ghosts()) return;

		auto map = vec_->getMap();
		auto ghost_map = ghosted_vec_->getMap();

		Tpetra::Export<
			LO,
			GO,
			vector_type::node_type> exporter(ghost_map, map);


		Teuchos::RCP<vector_type> y(new vector_type(map, 1));

		y->doExport(*ghosted_vec_, exporter, Tpetra::ADD);

		Tpetra::Import<vector_type::local_ordinal_type,
		                vector_type::global_ordinal_type,
		                vector_type::node_type> importer(map, ghost_map);

		 ghosted_vec_->doImport(*y, importer, Tpetra::INSERT);
	}

	TpetraVector::TpetraVector(const TpetraVector &other)
	{ 
		copy(other);
	}

	void TpetraVector::copy(const TpetraVector &other)
	{
		if(other.has_ghosts()) {
			ghosted_vec_ = Teuchos::rcp(new vector_type(other.ghosted_vec_->getMap(), 1));
			ghosted_vec_->assign(*other.ghosted_vec_);
			vec_ = ghosted_vec_->offsetViewNonConst(other.vec_->getMap(), 0);
		} else {
			if(vec_.is_null() || other.size().get(0) != size().get(0)) {
				vec_ = (Teuchos::rcp(new vector_type(*other.vec_, Teuchos::Copy)));
			} else {
				vec_->assign(other.implementation());
			}
		}
	}

	TpetraVector &TpetraVector::operator=(const TpetraVector &other)
	{
	    if(this == &other) return *this;

	    if(other.is_null()) {
	        vec_.reset();
	        return *this;
	    }

	    copy(other);
	    return *this;
	}


	void TpetraVector::write_unlock(WriteMode mode)
	{
	    switch(mode) {
	        case utopia::GLOBAL_ADD: {
	            export_ghosts_add();
	            break;
	        }

	        case utopia::LOCAL: {
	        	break;
	        }

	        default: {
	        	update_ghosts(); 
	            break;
	        }
	    }

	    write_data_ = Teuchos::ArrayRCP<Scalar>();
	}

}

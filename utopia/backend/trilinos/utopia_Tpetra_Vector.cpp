#include "utopia_Tpetra_Vector.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Instance.hpp"

#include <Tpetra_CrsMatrix_decl.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <cmath>

namespace utopia {

	void TpetraVector::add_vector(
	    const std::vector<global_ordinal_type> &indices,
	    const std::vector<Scalar> &values)
	{
		const std::size_t n = values.size();
		assert(n == indices.size());

		for(std::size_t i = 0; i < n; ++i) {
			add(indices[i], values[i]);
		}
	}

	void TpetraVector::set_vector(
	    const std::vector<global_ordinal_type> &indices,
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

	TpetraVector::Scalar TpetraVector::sum() const
	{
	    m_utopia_warning_once("> TpetraVector::sum is hand-coded");

	    auto data = implementation().getData();

	    Scalar ret_temp = 0.;

	    for(auto i = 0; i < data.size(); ++i) {
	        ret_temp += data[i];
	    }

	    double ret = ret_temp;
	    auto &comm = *communicator();
	    double ret_global = 0.;

	    Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 1, &ret, &ret_global);
	    return ret_global;
	}

	TpetraVector::Scalar TpetraVector::min() const
	{
	    m_utopia_warning_once("> TpetraVector::min is hand-coded");

	    auto data = implementation().getData();

	    Scalar ret_temp = data[0];

	    for(auto i = 1; i < data.size(); ++i) {
	        ret_temp = std::min(data[i], ret_temp);
	    }

	    double ret = ret_temp;
	    auto &comm = *communicator();
	    double ret_global = 0.;

	    Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, 1, &ret, &ret_global);
	    return ret_global;
	}

	TpetraVector::Scalar TpetraVector::max() const
	{
	    m_utopia_warning_once("> TpetraVector::max is hand-coded");

	    auto data = implementation().getData();

	    Scalar ret_temp = data[0];

	    for(auto i = 1; i < data.size(); ++i) {
	        ret_temp = std::max(data[i], ret_temp);
	    }

	    double ret = ret_temp;
	    auto &comm = *communicator();
	    double ret_global = 0.;

	    Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &ret, &ret_global);
	    return ret_global;
	}

	bool TpetraVector::is_nan_or_inf() const
	{
		m_utopia_warning_once("> TpetraVector::is_nan_or_inf is hand-coded");

		auto data = implementation().getData();

		int ret = 0;

		for(auto i = 0; i < data.size(); ++i) {
		    if(std::isnan(data[i]) || std::isinf(data[i])) {
		    	ret = 1;
		    	break;
		    }
		}

		auto &comm = *communicator();
		int ret_global = 0;

		Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &ret, &ret_global);
		return ret_global;
	}

	void TpetraVector::ghosted(
		const rcp_comm_type &comm, 
		const TpetraVector::global_ordinal_type &local_size,
	    const TpetraVector::global_ordinal_type &global_size,
		const std::vector<global_ordinal_type> &ghost_index
	)
	{
		rcp_map_type map(new map_type(global_size, local_size, 0, comm));
		rcp_map_type ghost_map;

		if(!ghost_index.empty()) {

			Range r = { map->getMinGlobalIndex(), map->getMaxGlobalIndex() + 1 };

			std::vector<global_ordinal_type> filled_with_local;
			filled_with_local.reserve(r.extent() + ghost_index.size());

			for(auto i = r.begin(); i != r.end(); ++i) {
				filled_with_local.push_back(i);
			}

			for(auto g : ghost_index) {
				if(!r.inside(g)) {
					filled_with_local.push_back(g);
				}
			}

			const Teuchos::ArrayView<const global_ordinal_type>
			   local_indices(filled_with_local);

			 ghost_map = Teuchos::rcp(new map_type(global_size, local_indices, 0, comm));

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
			local_ordinal_type,
			global_ordinal_type,
			vector_type::node_type> importer(map, ghost_map);

		ghosted_vec_->doImport(*vec_, importer, Tpetra::INSERT);
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
			if(vec_.is_null() == other.size().get(0) != size().get(0)) {
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
}

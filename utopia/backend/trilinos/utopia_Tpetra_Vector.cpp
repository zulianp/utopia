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
}

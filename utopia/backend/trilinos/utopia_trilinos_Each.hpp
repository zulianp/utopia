#ifndef UTOPIA_TRILINOS_EACH_HPP
#define UTOPIA_TRILINOS_EACH_HPP

#include "utopia_trilinos_Types.hpp"
#include "utopia_For.hpp"

namespace utopia {

	template<class Fun>
	inline void each_apply(TSMatrixd &mat, Fun fun) 
	{	
		auto r = row_range(mat);

		if(r.empty()) {
			return;
		}

		auto impl = raw_type(mat);
		auto local_mat = impl->getLocalMatrix();
		auto n = local_mat.numRows();

		for(decltype(n) i = 0; i < n; ++i) {
			auto row = local_mat.row(i);
			auto n_values = row.length;
			
			for(decltype(n_values) k = 0; k < n_values; ++k) {
				auto &val = row.value(k);
				val = fun(val);
			}
		}
	}

	template<class Fun>
	inline void each_transform(TSMatrixd &mat, Fun fun) 
	{	
		auto rr = row_range(mat);
		if(rr.empty()) return;

		auto impl = raw_type(mat);
		auto col_map   = impl->getColMap()->getLocalMap();
		auto row_map   = impl->getRowMap()->getLocalMap();
		auto local_mat = impl->getLocalMatrix();

		auto n = local_mat.numRows();

		for(decltype(n) i = 0; i < n; ++i) {
			auto row = local_mat.row(i);
			auto n_values = row.length;
			
			for(decltype(n_values) k = 0; k < n_values; ++k) {
				auto &val = row.value(k);
				const auto global_row = row_map.getGlobalElement(i);
				const auto global_col = col_map.getGlobalElement(row.colidx(k));
				
				val = fun(global_row, global_col, val);
			}
		}
	}

	template<class Fun>
	inline void each_read(const TVectord &v, Fun fun) 
	{
		auto view = raw_type(v)->getLocalView<Kokkos::HostSpace>();

		const auto r = range(v);
		const std::size_t r_begin = r.begin();
		
		For<>::apply(
			r_begin,
			r.end(),
			[r_begin, &fun, &view](const std::size_t i) {
				fun(i, view(i - r_begin, 0));
			}
		);
	}

	template<class Fun>
	inline void each_write(TVectord &v, Fun fun) 
	{
		auto view = raw_type(v)->getLocalView<Kokkos::HostSpace>();

		const auto r = range(v);
		const std::size_t r_begin = r.begin();
		
		For<>::apply(
			r_begin,
			r.end(),
			[r_begin, &fun, &view](const std::size_t i) {
				view(i - r_begin, 0) = fun(i);
			}
		);
	}

	template<class Fun>
	inline void each_transform(const TVectord &in, TVectord &out, Fun fun) 
	{
		const auto r = range(out);
		const std::size_t r_begin = r.begin();

		if(&in == &out) {
			auto view = raw_type(out)->getLocalView<Kokkos::HostSpace>();

			For<>::apply(
				r_begin,
				r.end(),
				[r_begin, &fun, &view](const std::size_t i) {
					const auto idx = i - r_begin;
					auto &val = view(idx, 0);
					val = fun(i, val);
				}
			);
		} else {

			auto view_in  = raw_type(in)->getLocalView<Kokkos::HostSpace>();
			auto view_out = raw_type(out)->getLocalView<Kokkos::HostSpace>();

			For<>::apply(
				r_begin,
				r.end(),
				[r_begin, &fun, &view_in, &view_out](const std::size_t i) {
					const auto idx = i - r_begin;
					const auto &val = view_in(idx, 0);
					view_out(idx, 0) = fun(i, val);
				}
			);
		}
	}

}

#endif //UTOPIA_TRILINOS_EACH_HPP

#ifndef UTOPIA_TRILINOS_EACH_HPP
#define UTOPIA_TRILINOS_EACH_HPP

#include "utopia_trilinos_Types.hpp"

namespace utopia {

	template<class Fun>
	inline void each_apply(TSMatrixd &mat, Fun fun) 
	{	
		auto r = row_range(mat);

		if(r.empty()) {
			return;
		}

		auto impl = raw_type(mat);
		auto col_map = impl->getColMap();
		auto row_map = impl->getRowMap();
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
		auto col_map = impl->getColMap();
		auto row_map = impl->getRowMap();
		auto local_mat = impl->getLocalMatrix();

		auto n = local_mat.numRows();

		for(decltype(n) i = 0; i < n; ++i) {
			auto row = local_mat.row(i);
			auto n_values = row.length;
			
			for(decltype(n_values) k = 0; k < n_values; ++k) {
				auto &val = row.value(k);
				const auto global_row = row_map->getGlobalElement(i);
				const auto global_col = col_map->getGlobalElement(row.colidx(k));
				
				val = fun(global_row, global_col, val);
			}
		}
	}
		
}

#endif //UTOPIA_TRILINOS_EACH_HPP

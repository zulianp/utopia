#ifndef UTOPIA_TRILINOS_EACH_PARALLEL_HPP
#define UTOPIA_TRILINOS_EACH_PARALLEL_HPP

#include "Kokkos_Core.hpp"
#include "utopia_trilinos_Types.hpp"

namespace utopia {
    typedef Kokkos::TeamPolicy<>               team_policy;
	typedef Kokkos::TeamPolicy<>::member_type  member_type;

	template<class Fun>
	inline void each_apply_parallel(TSMatrixd &mat, Fun fun) 
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


	    Kokkos::parallel_for( team_policy( n, Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type &teamMember) {
	    	
	    	const int j = teamMember.league_rank();
	    	auto row = local_mat.row(j);
	    	auto n_values = row.length;

	    	Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember, n_values), [&] (const int i) {
	    		auto &val = row.value(i);
	    		val = fun(val);
	    	});
	    });
	}
		
}

#endif //UTOPIA_TRILINOS_EACH_HPP

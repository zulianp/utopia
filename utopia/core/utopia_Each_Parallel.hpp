#ifndef UTOPIA_EACH_PARALLEL_HPP
#define UTOPIA_EACH_PARALLEL_HPP 

#include "utopia_ForwardDeclarations.hpp"

#include "utopia_Base.hpp"

#include "utopia_Each.hpp"


#include <Tpetra_Core.hpp>
// This is the only header file you need to include for the "core"
// part of Kokkos.  That includes Kokkos::View, Kokkos::parallel_*,
// and atomic updates.
#include <Kokkos_Core.hpp>

#include <memory>


namespace utopia {


	template<class Tensor, int Order = Tensor::Order, int FILL_TYPE = Tensor::FILL_TYPE>
	class ParallelEach {};

	template<int FILL_TYPE>
    class ParallelEach<TVectord, 1, FILL_TYPE>{
	public:
		template<class Fun>
		inline static void apply_write_parallel(const TVectord &v, Fun fun)
		{

            auto k_v = v.implementation().implementation().template getLocalView<Kokkos::HostSpace> ();
            Kokkos::parallel_for (k_v.extent(0), KOKKOS_LAMBDA (const int i) {
                k_v(i,0) = fun(i);
            });
        }

    template<class Fun>
		inline static void apply_read_parallel(const TVectord &v, Fun fun)
		{

            auto k_v = v.implementation().implementation().template getLocalView<Kokkos::HostSpace> ();
            Kokkos::parallel_for (k_v.extent(0), KOKKOS_LAMBDA (const int i) {
                fun(i,k_v(i,0));
            });
    }
		
};

        
    template<class Tensor, class Fun> 
	inline void each_read_parallel(const Tensor &v, Fun fun) 
	{
		ParallelEach<Tensor>::apply_read_parallel(v, fun);
	}

	template<class Tensor, class Fun>
	inline void each_write_parallel(Tensor &v, Fun fun) 
	{
		ParallelEach<Tensor>::apply_write_parallel(v, fun);
	}




}
#endif

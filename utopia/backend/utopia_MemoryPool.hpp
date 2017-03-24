#ifndef UTOPIA_MEMORY_POOL_HPP
#define UTOPIA_MEMORY_POOL_HPP

#include <map>

#include "utopia_Base.hpp"
#include "utopia_Size.hpp"

#ifdef WITH_PETSC
#include "petscmat.h"
#include "petscvec.h"
#endif // WITH_PETSC

namespace std {
	template<>
	struct less<std::pair<utopia::Size, utopia::Size>> {
		int operator()(const std::pair<utopia::Size, utopia::Size>& a, const std::pair<utopia::Size, utopia::Size>& b) const {
			if (a.second.n_dims() == b.second.n_dims()) {
				for (size_t i = 0; i < a.second.n_dims(); i++) {
					if (a.second.get(i) != b.second.get(i))
						return a.second.get(i) < b.second.get(i);
				}
				if (a.first.n_dims() == b.first.n_dims()) {
					for (size_t i = 0; i < a.first.n_dims(); i++) {
						if (a.first.get(i) != b.first.get(i))
							return a.first.get(i) < b.first.get(i);
					}
					return 0;
				} else {
					return a.first.n_dims() < b.first.n_dims();
				}
			} else {
				return a.second.n_dims() < b.second.n_dims();
			}
		}
	};
}

namespace utopia {

	class MemoryPool {
	public:
		static MemoryPool& getInstance();

		void fullGC();

#ifdef WITH_PETSC
		void setCommunicator(MPI_Comm);

		Vec* getVec(const Size& local, const Size& global);
		Vec* getVec(const Vec v);
		Mat* getMat(const Size& local, const Size& global);
		Mat* getMat(const Mat m);
		Mat* getSparseMat(const Size& local, const Size& global);
		Mat* getSparseMat(const Mat m);

		void put(Vec* v);
		void put(Mat* m);
		void putSparse(Mat* m);
#endif // WITH_PETSC

	private:
		MemoryPool() = default;
		~MemoryPool();

#ifdef WITH_PETSC
		MPI_Comm comm_;
		// TODO put a limit on maximum memory used
		std::multimap<std::pair<Size, Size>, Vec*> vec_pool_;
		std::multimap<std::pair<Size, Size>, Mat*> mat_pool_;
		std::multimap<std::pair<Size, Size>, Mat*> sparse_mat_pool_;
#endif // WITH_PETSC

	};


	#define MEMPOOL() MemoryPool::getInstance()
}

#endif // UTOPIA_MEMORY_POOL_HPP

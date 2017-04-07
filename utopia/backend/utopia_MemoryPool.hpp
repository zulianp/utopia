#ifndef UTOPIA_MEMORY_POOL_HPP
#define UTOPIA_MEMORY_POOL_HPP

#include <unordered_map>

#include "utopia_Base.hpp"
#include "utopia_Size.hpp"

#ifdef WITH_PETSC
#include "petscmat.h"
#include "petscvec.h"
#endif // WITH_PETSC

namespace std {
	template<>
	struct hash<std::pair<utopia::Size, utopia::Size>> {
		uint64_t operator()(const std::pair<utopia::Size, utopia::Size>& x) const {
			if (x.second.n_dims() == 0)
				return 0;
			if (x.second.n_dims() == 1)
				return (x.second.get(0) << 32) ^ (x.first.get(0));
			return ((x.second.get(0) ^ (x.second.get(1) << 12)) << 32)
				^ (x.first.get(0) ^ (x.first.get(1) << 12));
		}
	};

	template<>
	struct equal_to<std::pair<utopia::Size, utopia::Size>> {
		bool operator()(const std::pair<utopia::Size, utopia::Size>& a, const std::pair<utopia::Size, utopia::Size>& b) const {
			if (a.second.n_dims() == b.second.n_dims()) {
				for (size_t i = 0; i < a.second.n_dims(); i++) {
					if (a.second.get(i) != b.second.get(i))
						return false;
				}
				if (a.first.n_dims() == b.first.n_dims()) {
					for (size_t i = 0; i < a.first.n_dims(); i++) {
						if (a.first.get(i) != b.first.get(i))
							return false;
					}
					return true;
				} else {
					return false;
				}
			} else {
				return false;
			}
		}
	};
}

namespace utopia {

	class MemoryPool {
	public:
		static MemoryPool& getInstance();

		void GC();
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

		std::unordered_map<std::pair<Size, Size>, std::vector<Vec*>> vec_pool_;
		std::unordered_map<std::pair<Size, Size>, std::vector<Mat*>> mat_pool_;
		// std::unordered_map<std::pair<Size, Size>, std::vector<Mat*>> sparse_mat_pool_;
#endif // WITH_PETSC

	};


	#define MEMPOOL() MemoryPool::getInstance()
}

#endif // UTOPIA_MEMORY_POOL_HPP

#include "utopia_MemoryPool.hpp"

namespace utopia {

#ifdef WITH_PETSC
	static Size size(const Vec v) {
		PetscInt n;
		VecGetSize(v, &n);
		return {n};
	}

	static Size size(const Mat mat) {
		PetscInt m, n;
		MatGetSize(mat, &m, &n);
		return {m, n};
	}

	static Size local_size(const Vec v) {
		PetscInt n;
		VecGetLocalSize(v, &n);
		return {n};
	}

	static Size local_size(const Mat mat) {
		PetscInt m, n;
		MatGetLocalSize(mat, &m, &n);
		return {m, n};
	}
#endif //WITH_PETSC


	MemoryPool& MemoryPool::getInstance() {
		static MemoryPool instance;
		return instance;
	}

	void MemoryPool::fullGC() {
#ifdef WITH_PETSC
		for (auto it = vec_pool_.begin(); it != vec_pool_.end(); it++) {
			for (size_t i = 0; i < it->second.size(); i++) {
				VecDestroy(it->second[i]);
				delete it->second[i];
			}
		}
		vec_pool_.clear();
		for (auto it = mat_pool_.begin(); it != mat_pool_.end(); it++) {
			for (size_t i = 0; i < it->second.size(); i++) {
				MatDestroy(it->second[i]);
				delete it->second[i];
			}
		}
		mat_pool_.clear();
/*
		for (auto it = sparse_mat_pool_.begin(); it != sparse_mat_pool_.end(); it++) {
			MatDestroy(it->second);
			delete it->second;
		}
		sparse_mat_pool_.clear();
*/
#endif // WITH_PETSC
	}

	void MemoryPool::GC() {
#ifdef WITH_PETSC
		for (auto it = vec_pool_.begin(); it != vec_pool_.end(); it++) {
			if (it->second.size() == 0) {
				vec_pool_.erase(it);
			}
		}
		for (auto it = mat_pool_.begin(); it != mat_pool_.end(); it++) {
			if (it->second.size() == 0) {
				mat_pool_.erase(it);
			}
		}
/*
		for (auto it = mat_pool_.begin(); it != mat_pool_.end(); it++) {
			MatDestroy(it->second);
			delete it->second;
		}
		mat_pool_.clear();
		for (auto it = sparse_mat_pool_.begin(); it != sparse_mat_pool_.end(); it++) {
			MatDestroy(it->second);
			delete it->second;
		}
		sparse_mat_pool_.clear();
*/
#endif // WITH_PETSC
	}

#ifdef WITH_PETSC
	void MemoryPool::setCommunicator(MPI_Comm comm) {
		comm_ = comm;
	}

	Vec* MemoryPool::getVec(const Size& local, const Size& global) {
		// TODO consider PETSC_DECIDE and PETSC_DETERMINE in pool lookup
		auto it = vec_pool_.find(std::make_pair(local, global));
		if (it != vec_pool_.end() && it->second.size() > 0) {
			// std::cout << "Succesfully reused Vec" << std::endl;
			Vec* outval = *it->second.rbegin();
			it->second.pop_back();
			return outval;
		}

		// Allocation if not in pool
		// std::cout << "Had to allocate Vec" << std::endl;
		Vec* v = new Vec;
		if (global.n_dims() >= 1 && local.n_dims() >= 1) {
			VecCreateMPI(comm_, local.get(0), global.get(0), v);
		} else {
			VecCreate(comm_, v);
		}
		return v;
	}

	Vec* MemoryPool::getVec(const Vec v) {
		Size local = local_size(v), global = size(v);
		return getVec(local, global);
	}

	Mat* MemoryPool::getMat(const Size& local, const Size& global) {
		// TODO consider PETSC_DECIDE and PETSC_DETERMINE in pool lookup
		auto it = mat_pool_.find(std::make_pair(local, global));
		if (it != mat_pool_.end() && it->second.size() > 0) {
			// std::cout << "Succesfully reused Mat" << std::endl;
			Mat* outval = *it->second.rbegin();
			it->second.pop_back();
			return outval;
		}

		// Allocation if not in pool
		// std::cout << "Had to allocate Mat" << std::endl;
		Mat* m = new Mat;
		if (global.n_dims() >= 2 && local.n_dims() >= 2) {
			MatCreateDense(comm_, local.get(0), local.get(1), global.get(0), global.get(1), NULL, m);
		} else {
			MatCreate(comm_, m);
		}
		return m;
	}

	Mat* MemoryPool::getMat(const Mat v) {
		Size local = local_size(v), global = size(v);
		return getMat(local, global);
	}

	Mat* MemoryPool::getSparseMat(const Size& local, const Size& global) {
		// TODO consider PETSC_DECIDE and PETSC_DETERMINE in pool lookup
		/* Disabled for performance comparison
		auto it = sparse_mat_pool_.find(std::make_pair(local, global));
		if (it != sparse_mat_pool_.end()) {
			// std::cout << "Succesfully reused Sparse Mat" << std::endl;
			Mat* outval = it->second;
			sparse_mat_pool_.erase(it);
			return outval;
		} */

		// Allocation if not in pool
		// std::cout << "Had to allocate Sparse Mat" << std::endl;
		Mat* m = new Mat;
		if (global.n_dims() >= 2 && local.n_dims() >= 2) {
			MatCreateAIJ(comm_, local.get(0), local.get(1), global.get(0), global.get(1),
				1, PETSC_NULL, 1, PETSC_NULL, m);
		} else {
			MatCreate(comm_, m);
		}
		return m;
	}

	Mat* MemoryPool::getSparseMat(const Mat v) {
		Size local = local_size(v), global = size(v);
		return getSparseMat(local, global);
	}

	void MemoryPool::put(Vec* v) {
		PetscInt local, global;
		VecType t;
		if (VecGetType(*v, &t) || t == 0 || VecGetLocalSize(*v, &local) || VecGetSize(*v, &global)) {
			VecDestroy(v); // FIXME how can we reuse these?
			delete v;
		} else {
			auto k = std::make_pair<Size, Size>({local}, {global});
			if (vec_pool_[k].size() > 5) {
				VecDestroy(v);
				delete v;
			} else {
				vec_pool_[k].push_back(v);
			}
		}
	}

	void MemoryPool::put(Mat* m) {
		PetscInt local_m, local_n, global_m, global_n;
		MatType t;
		if (MatGetType(*m, &t) || t == 0 || MatGetLocalSize(*m, &local_m, &local_n) || MatGetSize(*m, &global_m, &global_n)) {
			MatDestroy(m); // FIXME how can we reuse these?
			delete m;
		} else {
			if (strcmp(t, "seqaij") == 0) {
				assert(false && "Disposing a sparse matrix like it's a dense matrix");
			}
			auto k = std::make_pair<Size, Size>({local_m, local_n}, {global_m, global_n});
			if (mat_pool_[k].size() > 5) {
				MatDestroy(m);
				delete m;
			} else {
				mat_pool_[k].push_back(m);
			}
		}
	}

	void MemoryPool::putSparse(Mat* m) {
		// PetscInt local_m, local_n, global_m, global_n;
		// MatType t;
		// if (MatGetType(*m, &t) || t == 0 || MatGetLocalSize(*m, &local_m, &local_n) || MatGetSize(*m, &global_m, &global_n)) {
			MatDestroy(m); // FIXME how can we reuse these?
			delete m;
		// } else {
		// 	sparse_mat_pool_.insert(std::make_pair(std::make_pair<Size, Size>({local_m, local_n}, {global_m, global_n}), m));
		// }
	}

#endif // WITH_PETSC

	MemoryPool::~MemoryPool() {
		// FIXME causes error "Call to MPI after MPI has been finalized"
		// gc();
	}
}


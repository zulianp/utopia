#include "utopia_MemoryPool.hpp"

std::less<std::pair<utopia::Size, utopia::Size>> comparator;


namespace utopia {

	MemoryPool& MemoryPool::getInstance() {
		static MemoryPool instance;
		return instance;
	}

#ifdef WITH_PETSC
	void MemoryPool::setCommunicator(MPI_Comm comm) {
		comm_ = comm;
	}

	Vec* MemoryPool::getVec(const Size& local, const Size& global) {
		auto it = vec_pool_.find(std::make_pair(local, global));
		if (it != vec_pool_.end()) {
			std::cout << "Succesfully reused Vec" << std::endl;
			Vec* outval = it->second;
			vec_pool_.erase(it);
			return outval;
		}

		// Allocation if not in pool
		Vec* v = new Vec;
		if (global.n_dims() >= 1 && local.n_dims() >= 1) {
			VecCreateMPI(comm_, local.get(0), global.get(0), v);
		} else {
			VecCreate(comm_, v);
		}
		return v;
	}

	Mat* MemoryPool::getMat(const Size& local, const Size& global) {
		auto it = mat_pool_.find(std::make_pair(local, global));
		if (it != mat_pool_.end()) {
			std::cout << "Succesfully reused Mat" << std::endl;
			Mat* outval = it->second;
			mat_pool_.erase(it);
			return outval;
		}

		// Allocation if not in pool
		Mat* m = new Mat;
		if (global.n_dims() >= 2 && local.n_dims() >= 2) {
			MatCreateDense(comm_, local.get(0), local.get(1), global.get(0), global.get(1), NULL, m);
		} else {
			MatCreate(comm_, m);
		}
		return m;
	}

	Mat* MemoryPool::getSparseMat(const Size& local, const Size& global) {
		auto it = sparse_mat_pool_.find(std::make_pair(local, global));
		if (it != sparse_mat_pool_.end()) {
			Mat* outval = it->second;
			sparse_mat_pool_.erase(it);
			return outval;
		}

		// Allocation if not in pool
		Mat* m = new Mat;
		if (global.n_dims() >= 2 && local.n_dims() >= 2) {
			MatCreateAIJ(comm_, local.get(0), local.get(1), global.get(0), global.get(1),
				1, PETSC_NULL, 1, PETSC_NULL, m);
		} else {
			MatCreate(comm_, m);
		}
		return m;
	}

	void MemoryPool::put(Vec* v, Size& local, Size& global) {
		vec_pool_.insert(std::make_pair(std::make_pair(local, global), v));
	}

	void MemoryPool::put(Mat* m, Size& local, Size& global) {
		mat_pool_.insert(std::make_pair(std::make_pair(local, global), m));
	}

	void MemoryPool::putSparse(Mat* m, Size& local, Size& global) {
		sparse_mat_pool_.insert(std::make_pair(std::make_pair(local, global), m));
	}

#endif // WITH_PETSC

	MemoryPool::MemoryPool() : vec_pool_(comparator), mat_pool_(comparator), sparse_mat_pool_(comparator) {}

	MemoryPool::~MemoryPool() {

#ifdef WITH_PETSC
		for (auto it = vec_pool_.begin(); vec_pool_.begin() != vec_pool_.end(); it++) {
			VecDestroy(it->second);
			delete it->second;
		}
		for (auto it = mat_pool_.begin(); mat_pool_.begin() != mat_pool_.end(); it++) {
			MatDestroy(it->second);
			delete it->second;
		}
		for (auto it = sparse_mat_pool_.begin(); sparse_mat_pool_.begin() != sparse_mat_pool_.end(); it++) {
			MatDestroy(it->second);
			delete it->second;
		}
#endif // WITH_PETSC

	}
}


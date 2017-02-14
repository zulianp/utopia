#ifndef UTOPIA_MEMORY_HPP
#define UTOPIA_MEMORY_HPP

#include <memory>
#include "petscmat.h"
#include "petscvec.h"

#include "utopia_Core.hpp"

namespace utopia {
	
	template <typename T, int FillType>
	class Allocator {
	public:
		static std::shared_ptr<T> claim(Size global, Size local) {
			assert(false && "No Allocator known for this type");
			return nullptr;
		}
	};

	template <typename T, int FillType>
	class Memory {
	public:
		Memory(Size global, Size local) : is_owner_(true) {
			mem_ = Allocator<T, FillType>::claim(global, local);
		}

		Memory(const Memory& m) : is_owner_(false), mem_(m.mem_) { }

		Memory(Memory&& m) : is_owner_(m.is_owner_), mem_(std::move(m.mem_)) {
			m.is_owner_ = false;
		}

		T& implementation() {
			return *mem_;
		}
		
	private:
		bool is_owner_;
		std::shared_ptr<T> mem_;
	};
	
	
	template<>
	class Allocator<Mat, 0> {
	public:
		static std::shared_ptr<Mat> claim(Size global, Size local) {
			assert(global.n_dims() == 2);
			assert(local.n_dims() == 2);
			
			//TODO ask pool for m
			Mat* m = new Mat;
			MatCreate(PETSC_COMM_WORLD, m);
			
			return std::shared_ptr<Mat>(m, [] (Mat* m){
				//TODO return m to pool
				MatDestroy(m);
				delete m;
			});
		}
	};

	template<>
	class Allocator<Vec, 0> {
	public:
		static std::shared_ptr<Vec> claim(Size global, Size local) {
			assert(global.n_dims() == 1);
			assert(local.n_dims() == 1);
			
			//TODO ask pool for v
			Vec* v = new Vec;
			VecCreateMPI(PETSC_COMM_WORLD, local.get(0), global.get(0), v);
			
			return std::shared_ptr<Vec>(v, [] (Vec* v){
				//TODO return v to pool
				VecDestroy(v);
				delete v;
			});
		}
	};

}

#endif //UTOPIA_MEMORY_HPP
// FIXME - only used with PETSc for now
#ifdef WITH_PETSC

#ifndef UTOPIA_MEMORY_HPP
#define UTOPIA_MEMORY_HPP

#include <memory>

#include "utopia_Base.hpp"
#include "utopia_Size.hpp"

namespace utopia {

	template <typename T>
	using MemoryPtr = std::unique_ptr<T, std::function<void(T*)>>;


	template <typename T, int FillType>
	class Allocator {
	public:
		static MemoryPtr<T> claim(MPI_Comm comm, Size local, Size global) {
			assert(false && "No Allocator known for this type");
			return nullptr;
		}

		static MemoryPtr<T> clone(const MemoryPtr<T>&) {
			assert(false && "No Allocator known for this type");
			return nullptr;
		}
	};


	template <typename T, int FillType>
	class Memory {
	public:
		Memory(MPI_Comm comm = PETSC_COMM_WORLD) : comm_(comm), init_(false), is_owner_(true) {
			mem_ = Allocator<T, FillType>::claim(comm, {}, {});
		}

		Memory(MPI_Comm comm, const Size& local, const Size& global) : comm_(comm), init_(true), is_owner_(true) {
			mem_ = Allocator<T, FillType>::claim(comm, local, global);
		}

		Memory(const Memory& m) : comm_(m.comm_), init_(m.init_), is_owner_(true) {
			mem_ = Allocator<T, FillType>::clone(m.mem_);
		}

		Memory(Memory&& m) : mem_(std::move(m.mem_)), comm_(m.comm_), init_(m.init_), is_owner_(true) { }

		Memory& operator=(const Memory& m) {
			if (!is_used_ && is_owner_)
				std::cerr << "[Warning] Destroying an unused object!" << std::endl;
			comm_ = m.comm_;
			init_ = m.init_;
			is_owner_ = true;
			is_used_ = false;
			mem_ = Allocator<T, FillType>::clone(m.mem_);
			return *this;
		}

		void wrap(T& t) {
			if (!is_used_ && is_owner_)
				std::cerr << "[Warning] Destroying an unused object!" << std::endl;
			mem_ = MemoryPtr<T>(&t, [](T*){});
			is_owner_ = false;
			is_used_ = false;
		}

		void resize(const Size& local, const Size& global) {
			if (!is_used_ && is_owner_)
				std::cerr << "[Warning] Destroying an unused object!" << std::endl;
			mem_ = Allocator<T, FillType>::claim(comm_, local, global);
			init_ = true;
			is_used_ = false;
		}

		T& implementation() {
			is_used_ = true;
			return *mem_;
		}

		const T& implementation() const {
			is_used_ = true;
			return *mem_;
		}

		MPI_Comm& communicator() {
			return comm_;
		}

		const MPI_Comm& communicator() const {
			return comm_;
		}

		//FIXME this doesn't actually change the communicator used by the underlying type
		void setCommunicator(const MPI_Comm comm) {
			comm_ = comm;
		}

		bool isInitialized() const {
			return init_;
		}

	private:
		MemoryPtr<T> mem_;
		MPI_Comm comm_;
		bool init_;
		bool is_owner_;
		mutable bool is_used_ = false;
	};

}

#endif //UTOPIA_MEMORY_HPP
#endif //WITH_PETSC
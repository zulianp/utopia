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
		Memory(MPI_Comm comm) : comm_(comm), init_(false) {
			mem_ = Allocator<T, FillType>::claim(comm, {}, {});
		}

		Memory(MPI_Comm comm, const Size& local, const Size& global) : comm_(comm), init_(true) {
			mem_ = Allocator<T, FillType>::claim(comm, local, global);
		}

		Memory(const Memory& m) : comm_(m.comm_), init_(m.init_) {
			mem_ = Allocator<T, FillType>::clone(m.mem_);
		}

		Memory(Memory&& m) : comm_(m.comm_), mem_(std::move(m.mem_)), init_(m.init_) { }

		Memory& operator=(const Memory& m) {
			comm_ = m.comm_;
			init_ = m.init_;
			mem_ = Allocator<T, FillType>::clone(m.mem_);
			return *this;
		}

		void initialize(const Size& local, const Size& global) {
			mem_ = Allocator<T, FillType>::claim(comm_, local, global);
			init_ = true;
		}

		void resize(const Size& local, const Size& global) {
			mem_ = Allocator<T, FillType>::claim(comm_, local, global);
			init_ = true;
		}

		void destroy() {
			mem_ = nullptr;
			init_ = false;
		}

		T& implementation() {
			return *mem_;
		}

		MPI_Comm& communicator() {
			return comm_;
		}

		const MPI_Comm& communicator() const {
			return comm_;
		}

		void setCommunicator(const MPI_Comm comm) {
			comm_ = comm;
		}

		const T& implementation() const {
			return *mem_;
		}

		bool isInitialized() const {
			return init_;
		}

	private:
		MemoryPtr<T> mem_;
		MPI_Comm comm_;
		bool init_;
	};

}

#endif //UTOPIA_MEMORY_HPP
#endif //WITH_PETSC
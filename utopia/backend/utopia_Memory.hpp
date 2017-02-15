// FIXME - only used with PETSc for now
#ifdef WITH_PETSC

#ifndef UTOPIA_MEMORY_HPP
#define UTOPIA_MEMORY_HPP

#include <memory>

#include "utopia_Base.hpp"
#include "utopia_Size.hpp"

namespace utopia {

	template <typename T, int FillType>
	class Allocator {
	public:
		static std::unique_ptr<T, std::function<void(T*)>> claim(MPI_Comm comm, Size local, Size global) {
			assert(false && "No Allocator known for this type");
			return nullptr;
		}

		static std::unique_ptr<T, std::function<void(T*)>> clone(const std::unique_ptr<T, std::function<void(T*)>>&) {
			assert(false && "No Allocator known for this type");
			return nullptr;
		}
	};

	template <typename T, int FillType>
	class Memory {
	public:
		Memory(MPI_Comm comm) : comm_(comm) {
			mem_ = nullptr;
		}
		
		Memory(MPI_Comm comm, Size local, Size global) : comm_(comm) {
			mem_ = Allocator<T, FillType>::claim(comm, local, global);
		}

		Memory(const Memory& m) : comm_(m.comm_) {
			mem_ = Allocator<T, FillType>::clone(m.mem_);
		}

		Memory(Memory&& m) : comm_(m.comm_), mem_(std::move(m.mem_)) { }
		
		Memory& operator=(const Memory& m) {
			comm_ = m.comm_;
			mem_ = Allocator<T, FillType>::clone(m.mem_);
			return *this;
		}
		
		void initialize(Size local, Size global) {
			mem_ = Allocator<T, FillType>::claim(comm_, local, global);
		}
		
		void resize(Size local, Size global) {
			mem_ = Allocator<T, FillType>::claim(comm_, local, global);
		}

		void destroy() {
			mem_ = nullptr;
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
			return static_cast<bool>(mem_);
		}

	protected:
		std::unique_ptr<T, std::function<void(T*)>> mem_; //TODO why were we using shared_ptr?
		MPI_Comm comm_;
	};

}

#endif //UTOPIA_MEMORY_HPP
#endif //WITH_PETSC
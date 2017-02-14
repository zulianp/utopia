#ifndef UTOPIA_MEMORY_HPP
#define UTOPIA_MEMORY_HPP

#include <memory>

#include "utopia_Base.hpp"
#include "utopia_Size.hpp"

namespace utopia {

	template <typename T, int FillType>
	class Allocator {
	public:
		static std::unique_ptr<T> claim(MPI_Comm comm, Size global, Size local) {
			assert(false && "No Allocator known for this type");
			return nullptr;
		}

		static std::unique_ptr<T> clone(std::unique_ptr<T>&) {
			assert(false && "No Allocator known for this type");
			return nullptr;
		}
	};

	template <typename T, int FillType>
	class Memory {
	public:
		Memory(MPI_Comm comm, Size global, Size local) : comm_(comm) {
			mem_ = Allocator<T, FillType>::claim(comm, global, local);
		}

		Memory(const Memory& m) : comm_(m.comm_) {
			mem_ = Allocator<T, FillType>::clone(m.mem_);
		}

		Memory(Memory&& m) : comm_(m.comm_), mem_(std::move(m.mem_)) { }
		
		Memory& operator=(const Memory& m) {
			comm_ = m.comm_;
			mem_ = Allocator<T, FillType>::clone(m.mem_);
		}

		T& implementation() {
			return *mem_;
		}

		MPI_Comm& commmunicator() {
			return comm_;
		}

		const MPI_Comm& commmunicator() const {
			return comm_;
		}

		void setCommunicator(const MPI_Comm comm) {
            comm_ = comm;
        }

		const T& implementation() const {
			return *mem_;
		}

		bool isInitialized() const { //legacy signature
			return static_cast<bool>(mem_);
		}

		bool is_initialized() const {
			return static_cast<bool>(mem_);
		}

	private:
		std::unique_ptr<T> mem_; //TODO why were we using shared_ptr?
		MPI_Comm comm_;
	};

}

#endif //UTOPIA_MEMORY_HPP
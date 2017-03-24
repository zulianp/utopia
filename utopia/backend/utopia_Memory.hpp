// FIXME - only used with PETSc for now
#ifdef WITH_PETSC

#ifndef UTOPIA_MEMORY_HPP
#define UTOPIA_MEMORY_HPP

#include <memory>

#include "utopia_MemoryPool.hpp"

namespace utopia {

	template <typename T>
	using MemoryPtr = std::unique_ptr<T, std::function<void(T*)>>;


	template <typename T, int FillType>
	class Allocator {
	public:
		static void destructor(T*) {
			assert(false && "No Allocator known for this type");
		}

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
		template<typename T2, int FillType2>
		friend class Memory;

		Memory(MPI_Comm comm = PETSC_COMM_WORLD) : mem_(nullptr), comm_(comm), init_(false), is_owner_(true) { }

		Memory(MPI_Comm comm, const Size& local, const Size& global) : comm_(comm), init_(true), is_owner_(true) {
			mem_ = Allocator<T, FillType>::claim(comm, local, global);
		}

		Memory(MPI_Comm comm, T& t, bool own) : comm_(comm), init_(true), is_owner_(own) {
			mem_ = MemoryPtr<T>(&t, own ? Allocator<T, FillType>::destructor : [](T*){});
		}

		Memory(const Memory& m) : comm_(m.comm_), init_(m.init_), is_owner_(true) {
			mem_ = Allocator<T, FillType>::clone(m.mem_);
		}

		Memory(Memory&& m) : mem_(std::move(m.mem_)), comm_(m.comm_), init_(m.init_), is_owner_(m.is_owner_) {
			m.init_ = false;
			m.is_owner_ = false;
		}

		Memory& operator=(const Memory& m) {
			if (!is_used_ && is_owner_ && mem_)
				std::cout << "[Warning] Destroying an unused object!" << std::endl;
			mem_ = Allocator<T, FillType>::clone(m.mem_);
			comm_ = m.comm_;
			init_ = m.init_;
			is_owner_ = true;
			is_used_ = false;
			return *this;
		}

		template<int FillTypeOther>
		Memory& operator=(const Memory<T, FillTypeOther>& m) {
			if (!is_used_ && is_owner_ && mem_)
				std::cout << "[Warning] Destroying an unused object!" << std::endl;
			mem_ = Allocator<T, FillType>::clone(m.mem_);
			comm_ = m.comm_;
			init_ = m.init_;
			is_owner_ = true;
			is_used_ = false;
			return *this;
		}

		void init(const Size& local, const Size& global) {
			if (!is_used_ && is_owner_ && mem_)
				std::cout << "[Warning] Destroying an unused object!" << std::endl;
			mem_ = Allocator<T, FillType>::claim(comm_, local, global);
			init_ = true;
			is_owner_ = true;
			is_used_ = false;
		}

		void init(T& t, bool own) {
			if (!is_used_ && is_owner_ && mem_)
				std::cout << "[Warning] Destroying an unused object!" << std::endl;
			mem_ = MemoryPtr<T>(&t, own ? Allocator<T, FillType>::destructor : [](T*){});
			init_ = true;
			is_owner_ = own;
			is_used_ = false;
		}

		void init() {
			init({}, {});
			init_ = false;
		}

		void wrap(T& t) {
			init(t, false);
		}

		void resize(const Size& local, const Size& global) {
			init(local, global);
		}

		T& implementation() {
			if (!mem_)
				init();
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
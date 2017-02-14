#ifndef UTOPIA_MEMORY_HPP
#define UTOPIA_MEMORY_HPP

#include <memory>

#include "utopia_Base.hpp"
#include "utopia_Size.hpp"

namespace utopia {

	template <typename T, int FillType>
	class Allocator {
	public:
		static std::unique_ptr<T> claim(Size global, Size local) {
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
		Memory(Size global, Size local) : is_owner_(true) {
			mem_ = Allocator<T, FillType>::claim(global, local);
		}

		Memory(const Memory& m) : is_owner_(false) {
			mem_ = Allocator<T, FillType>::clone(m.mem_);
		}

		Memory(Memory&& m) : is_owner_(m.is_owner_), mem_(std::move(m.mem_)) {
			m.is_owner_ = false;
		}

		T& implementation() {
			return *mem_;
		}

	private:
		bool is_owner_;
		std::unique_ptr<T> mem_; //TODO why were we using shared_ptr?
	};

}

#endif //UTOPIA_MEMORY_HPP
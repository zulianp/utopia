#ifndef UTOPIA_MEMORY_HPP
#define UTOPIA_MEMORY_HPP

#include <memory>

namespace utopia {
	
	// In specializations, use Traits to determine allocation type?
	template <typename T, int FillType>
	class Allocator {
	public:
		static std::shared_ptr<T> claim(Size global, Size local) {
			assert(false && "No Allocator known for this type");
			return nullptr;
		}

		static void release(std::shared_ptr<T> &mem) {
			assert(false && "No Allocator known for this type");
		}
	}

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

	private:
		bool is_owner_;
		std::shared_ptr<T> mem_;
	}

}

#endif //UTOPIA_MEMORY_HPP
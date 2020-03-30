#ifndef UTOPIA_LAYOUT_HPP
#define UTOPIA_LAYOUT_HPP

#include <utility>
#include <algorithm>

namespace utopia {

	template<class Comm, typename SizeType_, int Order>
	class Layout {
	public:
		using SizeType = SizeType_;

		inline SizeType &local_size(const SizeType i = 0)
		{
			assert(i < Order);
			return local_size_[i];
		}

		inline const SizeType &local_size(const SizeType i = 0) const
		{
			assert(i < Order);
			return local_size_[i];
		}

		inline SizeType &global_size(const SizeType i = 0)
		{
			assert(i < Order);
			return global_size_[i];
		}

		inline const SizeType &global_size(const SizeType i = 0) const
		{
			assert(i < Order);
			return global_size_[i];
		}

		const Comm &comm() const { return comm_; }

		template<typename...Args>
		Layout(const Comm &comm, Args &&... args)
		: comm_(comm)
		{
			init(std::forward<Args>(args)...);
		}

		template<typename OtherSizeType>
		Layout(const Layout<Comm, OtherSizeType, Order> &other)
		{
			std::copy(&other.local_size(0),  &other.local_size(0) + Order,  local_size_);
			std::copy(&other.global_size(0), &other.global_size(0) + Order, global_size_);
		}

		inline void init(const SizeType local_size[Order], SizeType global_size[Order])
		{
			std::copy(local_size,  local_size + Order,  local_size_);
			std::copy(global_size, global_size + Order, global_size_);
		}

		inline void init(const SizeType local_size, SizeType global_size)
		{
			local_size_[0] = local_size;
			global_size_[0] = global_size;
		}

	private:
		Comm comm_;
		SizeType local_size_[Order];
		SizeType global_size_[Order];
	};


	template<typename SizeType, class Comm>
	inline Layout<Comm, SizeType, 1> layout(const Comm &comm, const SizeType &local_size, const SizeType &global_size)
	{
		return Layout<Comm, SizeType, 1>(comm, local_size, global_size);
	}

	template<typename SizeType, class Comm>
	inline Layout<Comm, SizeType, 2> layout(const Comm &comm, const std::pair<SizeType, SizeType> &local_size, const std::pair<SizeType, SizeType> &global_size)
	{
		SizeType ls[2] = { local_size.first, local_size.second   };
		SizeType gs[2] = { global_size.first, global_size.second };
		return Layout<Comm, SizeType, 2>(comm, ls, gs);
	}

}

#endif //UTOPIA_LAYOUT_HPP

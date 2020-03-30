#ifndef UTOPIA_LAYOUT_HPP
#define UTOPIA_LAYOUT_HPP

#include <utility>
#include <algorithm>

#include "utopia_Traits.hpp"
#include "utopia_ForwardDeclarations.hpp"

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

		inline SizeType &size(const SizeType i = 0)
		{
			assert(i < Order);
			return size_[i];
		}

		inline const SizeType &size(const SizeType i = 0) const
		{
			assert(i < Order);
			return size_[i];
		}

		inline bool same_local_size(const Layout &other) const
		{
			for(int i = 0; i < Order; ++i) {
				if(local_size_[i] != other.local_size(i)) {
					return false;
				}
			}

			return true;
		}

		inline bool same_size(const Layout &other) const
		{
			for(int i = 0; i < Order; ++i) {
				if(size_[i] != other.size(i)) {
					return false;
				}
			}

			return true;
		}

		/// collective (communicator of this is used for the computation)
		inline bool same(const Layout &other) const
		{
			if(!comm().same(other.comm())) return false;
			if(!same_size(other)) return false;

			return comm().conjunction(same_local_size(other));
		}

		const Comm &comm() const { return comm_; }

		template<typename...Args>
		Layout(const Comm &comm, Args &&... args)
		: comm_(comm)
		{
			init(std::forward<Args>(args)...);
		}

		Layout()
		{
			for(int i = 0; i < Order; ++i) {
				local_size_[i] = 0;
				size_[i] = 0;
			}
		}

		template<typename OtherSizeType>
		Layout(const Layout<Comm, OtherSizeType, Order> &other)
		{
			std::copy(&other.local_size(0),  &other.local_size(0) + Order,  local_size_);
			std::copy(&other.size(0), &other.size(0) + Order, size_);
		}

		inline void init(const SizeType local_size[Order], SizeType size[Order])
		{
			std::copy(local_size,  local_size + Order,  local_size_);
			std::copy(size, size + Order, size_);
		}

		inline void init(const SizeType local_size, SizeType size)
		{
			local_size_[0] = local_size;
			size_[0] = size;
		}

	private:
		Comm comm_;
		SizeType local_size_[Order];
		SizeType size_[Order];
	};


	template<typename SizeType, class Comm>
	inline Layout<Comm, SizeType, 1> layout(const Comm &comm, const SizeType &local_size, const SizeType &size)
	{
		return Layout<Comm, SizeType, 1>(comm, local_size, size);
	}

	template<typename SizeType, class Comm>
	inline Layout<Comm, SizeType, 2> layout(const Comm &comm, const std::pair<SizeType, SizeType> &local_size, const std::pair<SizeType, SizeType> &size)
	{
		SizeType ls[2] = { local_size.first, local_size.second   };
		SizeType gs[2] = { size.first, size.second };
		return Layout<Comm, SizeType, 2>(comm, ls, gs);
	}

	template<class V>
	Layout<typename Traits<V>::Communicator, typename Traits<V>::SizeType, 1> layout(const Tensor<V, 1> &t)
	{
		const auto &vec = t.derived();
		return Layout<typename Traits<V>::Communicator, typename Traits<V>::SizeType, 1>(vec.comm(), vec.local_size(), vec.size());
	}

	// template<class V>
	// Layout<typename Traits<V>::Communicator, typename Traits<V>::SizeType, 1> row_layout(const Tensor<V, 2> &t)
	// {
	// 	const auto &mat = t.derived();
	// 	return Layout<typename Traits<V>::Communicator, typename Traits<V>::SizeType, 1>(mat.comm(), mat.local_size().get(0), mat.size().get(0));
	// }

	template<class V>
	Layout<typename Traits<V>::Communicator, typename Traits<V>::SizeType, 1> row_layout(const Operator<V> &op)
	{
		return Layout<typename Traits<V>::Communicator, typename Traits<V>::SizeType, 1>(op.comm(), op.local_size().get(0), op.size().get(0));
	}

}

#endif //UTOPIA_LAYOUT_HPP

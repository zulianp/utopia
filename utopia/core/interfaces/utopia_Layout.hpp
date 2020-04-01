#ifndef UTOPIA_LAYOUT_HPP
#define UTOPIA_LAYOUT_HPP

#include <utility>
#include <algorithm>
#include <cassert>

#include "utopia_Traits.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Communicator.hpp"
#include "utopia_Size.hpp"

namespace utopia {

	template<class Comm, int Order, typename LocalSizeType_, typename SizeType_ = LocalSizeType_ >
	class Layout {
	public:
		using SizeType 		= SizeType_;
		using LocalSizeType = LocalSizeType_;

		inline LocalSizeType_ &local_size(const int i = 0)
		{
			assert(i < Order);
			return local_size_[i];
		}

		inline const LocalSizeType_ &local_size(const int i = 0) const
		{
			assert(i < Order);
			return local_size_[i];
		}

		inline SizeType &size(const int i = 0)
		{
			assert(i < Order);
			return size_[i];
		}

		inline const SizeType &size(const int i = 0) const
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

		template<class OtherComm, typename OtherSizeType, typename OtherLocalSizeType>
		Layout(const Layout<OtherComm, Order, OtherLocalSizeType, OtherSizeType> &other)
		: comm_(other.comm())
		{
			std::copy(&other.local_size(0),  &other.local_size(0) + Order,  local_size_);
			std::copy(&other.size(0), &other.size(0) + Order, size_);
		}

		inline void init(const LocalSizeType local_size[Order], SizeType size[Order])
		{
			std::copy(local_size,  local_size + Order,  local_size_);
			std::copy(size, size + Order, size_);
		}

		inline void init(const Size &local, const Size &global)
		{
			auto &local_data  = local.data();
			auto &global_data = global.data();

			assert(int(local_data.size()) <= Order);
			assert(int(global_data.size()) <= Order);

			std::copy(local_data.begin(),  local_data.end(),    local_size_);
			std::copy(global_data.begin(),  global_data.end(),  size_);
		}

		inline void init(const LocalSizeType local_size, SizeType size)
		{
			local_size_[0] = local_size;
			size_[0] = size;
		}

	private:
		Comm comm_;
		LocalSizeType local_size_[Order];
		SizeType size_[Order];
	};

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//// SERIAL LAYOUTS //////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	template<typename SizeType>
	inline Layout<SelfCommunicator, 1, SizeType> serial_layout(const SizeType &size)
	{
		return Layout<SelfCommunicator, 1, SizeType>(SelfCommunicator(), size, size);
	}

	template<typename SizeType>
	inline Layout<SelfCommunicator, 2, SizeType> serial_layout(const SizeType &rows, const SizeType &cols)
	{
        SizeType ls[2] = { rows, cols  };
        SizeType gs[2]      = { rows, cols };
		return Layout<SelfCommunicator, 2, SizeType>(SelfCommunicator(), ls, gs);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//// PARALLEL LAYOUTS ////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	template<typename LocalSizeType, typename SizeType, class Comm>
	inline Layout<Comm, 1, LocalSizeType, SizeType> layout(const Comm &comm, const LocalSizeType &local_size, const SizeType &size)
	{
		return Layout<Comm, 1, LocalSizeType, SizeType>(comm, local_size, size);
	}

	template<typename LocalSizeType, typename SizeType, class Comm>
	inline Layout<Comm, 2, LocalSizeType, SizeType> layout(
		const Comm &comm,
		const LocalSizeType &local_rows,
		const LocalSizeType &local_cols,
		const SizeType &rows,
		const SizeType &cols)
	{
		LocalSizeType ls[2] = { local_rows, local_cols  };
		SizeType gs[2]	    = { rows, cols };

		return Layout<Comm, 2, LocalSizeType, SizeType>(comm, ls, gs);
	}

	template<class V>
	Layout<typename Traits<V>::Communicator, 1, typename Traits<V>::LocalSizeType, typename Traits<V>::SizeType> layout(const Tensor<V, 1> &t)
	{
		const auto &vec = t.derived();
		return Layout<typename Traits<V>::Communicator, 1, typename Traits<V>::LocalSizeType, typename Traits<V>::SizeType>(vec.comm(), vec.local_size(), vec.size());
	}

	template<class V>
	Layout<typename Traits<V>::Communicator, 1, typename Traits<V>::LocalSizeType, typename Traits<V>::SizeType> row_layout(const Operator<V> &op)
	{
		return Layout<typename Traits<V>::Communicator, 1, typename Traits<V>::LocalSizeType, typename Traits<V>::SizeType>(op.comm(), op.local_size().get(0), op.size().get(0));
	}

	template<class V>
	Layout<typename Traits<V>::Communicator, 1, typename Traits<V>::LocalSizeType, typename Traits<V>::SizeType> col_layout(const Operator<V> &op)
	{
		return Layout<typename Traits<V>::Communicator, 1, typename Traits<V>::LocalSizeType, typename Traits<V>::SizeType>(op.comm(), op.local_size().get(1), op.size().get(1));
	}


	template<class V>
	Layout<typename Traits<V>::Communicator, 2, typename Traits<V>::LocalSizeType, typename Traits<V>::SizeType> layout(const Operator<V> &op)
	{
		return Layout<typename Traits<V>::Communicator, 2, typename Traits<V>::LocalSizeType, typename Traits<V>::SizeType>(
			op.comm(),
			op.local_size(),
			op.size()
		);
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//// UTITLITY ///////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	template<typename LocalSizeType, typename SizeType, class Comm>
	inline Layout<Comm, 2, LocalSizeType, SizeType> square_matrix_layout(const Layout<Comm, 1, LocalSizeType, SizeType> &in)
	{
		LocalSizeType ls[2] = { in.local_size(), in.local_size()  };
		SizeType gs[2]	    = { in.size(), in.size() };
		return Layout<Comm, 2, LocalSizeType, SizeType>(in.comm(), ls, gs);
	}

}

#endif //UTOPIA_LAYOUT_HPP

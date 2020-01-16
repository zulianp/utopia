#ifndef UTOPIA_LAYOUT_HPP
#define UTOPIA_LAYOUT_HPP

namespace utopia {

	template<typename SizeType_, int Order>
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

	private:
		SizeType local_size_[Order];
		SizeType global_size_[Order];
	};

}

#endif //UTOPIA_LAYOUT_HPP

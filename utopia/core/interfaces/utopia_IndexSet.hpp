#ifndef UTOPIA_INDEX_SET_HPP
#define UTOPIA_INDEX_SET_HPP


#include "utopia_Enums.hpp"
#include "utopia_Range.hpp"
#include "utopia_Layout.hpp"

namespace utopia {

	template<typename Index_, typename SizeType_>
	class IndexSet {
	public:
		using Index    = Index_;
		using SizeType = SizeType_;

		virtual ~IndexSet() {}

		//locks
		virtual void read_lock()  = 0;
		virtual void write_lock() = 0;

		virtual void read_unlock()  = 0;
		virtual void write_unlock() = 0;

		//basic mutators
		virtual void set(const SizeType &i, const Index &value) = 0;
		virtual Index get(const SizeType &i) = 0;

		//print function
		virtual void describe() const = 0;

		//utility functions
		virtual bool empty() const = 0;
		virtual bool clear() const = 0;
	};

}


#endif //UTOPIA_INDEX_SET_HPP

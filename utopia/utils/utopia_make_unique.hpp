#ifndef UTOPIA_MAKE_UNIQUE_HPP
#define UTOPIA_MAKE_UNIQUE_HPP

#include <utility>
#include <memory>

namespace utopia {
	
	// note: this implementation does not disable this overload for array types
	template<typename T, typename... Args>
	std::unique_ptr<T> make_unique(Args&&... args)
	{
	    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
	}
}

#endif //UTOPIA_MAKE_UNIQUE_HPP

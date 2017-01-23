#ifndef UTOPIA_UTILS_HPP
#define UTOPIA_UTILS_HPP

#include <string>
#include <vector>
#include <memory>

#include "utopia_MPI.hpp"
#include "utopia_Chrono.hpp"

namespace utopia 
{



	bool read(const std::string &path, std::string &str);
	bool write(const std::string &path, const std::string &str);
	void disp(const std::string &str);
	void disp(const char * str);
	std::vector<std::string> split_string(const std::string &, char);

	/** \addtogroup Utilities
	 * @ingroup ptr_utilities
	 * @brief Provide basic functionalities for pointer manipulations.
	 *  @{
	 */

	template <typename T>
	struct EmptyDeleter
	{
	    EmptyDeleter() /* noexcept */
	    {
	    }

	    template <typename U>
	    EmptyDeleter(const EmptyDeleter<U>&,
	        typename std::enable_if<
	            std::is_convertible<U*, T*>::value
	        >::type* = nullptr) /* noexcept */
	    {
	    }

	    void operator()(T* const) const /* noexcept */
	    {
	        // do nothing
	    }
	};


	/*!
	 * @fn std::shared_ptr<T> make_ref(T &obj)
	 * @brief use only when you are sure that the scope of the pointer is limited to the life of obj
	 * @return a shared_ptr that does not delete the pointer 
	 */
	template<typename T>
	std::shared_ptr<T> make_ref(T &obj)
	{
	    return std::shared_ptr<T>(&obj, EmptyDeleter<T>());
	}

	/*!
	 * @fn std::shared_ptr<const T> make_ref(const T &obj)
	 * @brief use only when you are sure that the scope of the pointer is limited to the life of obj
	 * @return a shared_ptr that does not delete the pointer 
	 */
	template<typename T>
	std::shared_ptr<const T> make_ref(const T &obj)
	{
	    return std::shared_ptr<const T>(&obj, EmptyDeleter<const T>());
	}

	  /** @}*/
}

#endif //UTOPIA_UTILS_HPP

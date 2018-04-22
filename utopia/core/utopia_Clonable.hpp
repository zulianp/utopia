#ifndef UTOPIA_CLONABLE_HPP
#define UTOPIA_CLONABLE_HPP 

#include "utopia_Logger.hpp"

namespace utopia {
	class Clonable {
	public:
		virtual ~Clonable()
		{
			m_utopia_warning_once("Clone overrides are incomplete.");
		}

		/** @brief This method copies the relevant 
		 * settings but does not have to copy all the state variables
		 * such as buffers. Maybe we should change its name to somthing 
		 * more suitable ...
		 */
		virtual Clonable * clone() const = 0;

	};
}

#endif //UTOPIA_CLONABLE_HPP

#ifndef UTOPIA_KSP_TYPES_HPP
#define UTOPIA_KSP_TYPES_HPP


#include "utopia_Base.hpp"
#include "utopia_Utils.hpp"

#include <string>
#include <vector>


namespace utopia {

	class KSPTypes final {
	public:
	    inline static KSPTypes &instance()
	    {
	        static KSPTypes instance_;
	        return instance_;
	    }

	    inline const std::string &ksp(const SizeType i)
	    {
	        return ksp_.at(i);
	    }

	    inline const std::string &pc(const SizeType i)
	    {
	        return pc_.at(i);
	    }

	    inline const std::string &package(const SizeType i)
	    {
	        return package_.at(i);
	    }

	    inline bool is_ksp_valid(const std::string &type) const
	    {
	        return in_array(type, ksp_);
	    }

	    inline bool is_pc_valid(const std::string &type) const
	    {
	        return in_array(type, pc_);
	    }

	    inline bool is_solver_package_valid(const std::string &type) const
	    {
	        return in_array(type, package_);
	    }

	private:
	    KSPTypes();
	    std::vector<std::string> ksp_;              /*!< Valid options for direct solver types. */
	    std::vector<std::string> pc_;
	    std::vector<std::string> package_;
	};

}

#endif //UTOPIA_KSP_TYPES_HPP
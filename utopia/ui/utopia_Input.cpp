#include "utopia_Input.hpp"
#include <exception>
#include "utopia_ui.hpp"

#include "utopia_Instance.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Reporter.hpp"

namespace utopia {

    void Input::KeyDeprecated::key_exists() const {
        if (mpi_world_rank() == 0) {
            utopia::err() << "Input[Warning]: \"" << deprecated_key_
                          << "\" is deprecated and its usage will be removed in the future. Use \"" << new_key_
                          << "\" instead!\n ";
        }
    }

    void Input::KeyRequired::key_does_not_exists() const {
        if (mpi_world_rank() == 0) {
            utopia::err() << "Input[Error]: The input \"" << key_ << "\" is required! aborting program!\n";
        }
        assert(false);
        Utopia::Abort();
    }

    bool Configurable::import(const Path &path) {
        try {
            auto istr = open_istream(path.to_string());
            if (!istr) {
                return false;
            }

            read(*istr);
            return true;
        } catch (const std::exception &ex) {
            std::cerr << "[Error] " << ex.what() << std::endl;
            return false;
        }
    }

    bool Configurable::import(const std::string &key, const Path &path) {
        try {
            auto istr = open_istream(path.to_string());
            if (!istr) {
                return false;
            }

            istr->get(key, *this);
            return true;
        } catch (const std::exception &ex) {
            std::cerr << "[Error] " << ex.what() << std::endl;
            return false;
        }
    }

    void Configurable::print_usage(std::ostream &) const {}
    void Configurable::print_param_usage(std::ostream &os,
                                         const std::string &name,
                                         const std::string &type,
                                         const std::string &description,
                                         const std::string &default_settings) const {
        os << name << std::setw(25 - name.size()) << " : <" << type << ">" << std::setw(10 - type.size()) << std::right
           << " | " << description << std::setw(50 - description.size()) << std::right << " | " << default_settings
           << " \n";
    }
}  // namespace utopia

//
// Created by Patrick Zulian on 01/10/15.
//

#ifndef UTOPIA_UTOPIA_BACKENDINFO_HPP
#define UTOPIA_UTOPIA_BACKENDINFO_HPP

#include <string>
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {
    class BackendInfo {
    public:
        const std::string &get_name() const {
            return name_;
        }

        void set_name(const std::string &name) {
            this->name_ = name;
        }

        BackendInfo(const std::string &name = "undefined")
        : name_(name)
        {}

    private:
        std::string name_;
    };


    template<class T>
    auto backend_info(const T &)
    {
        return Traits<T>::backend_info();
    }
}

#endif //UTOPIA_UTOPIA_BACKENDINFO_HPP

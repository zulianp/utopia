//
// Created by Patrick Zulian on 01/10/15.
//

#ifndef UTOPIA_UTOPIA_BACKENDINFO_HPP
#define UTOPIA_UTOPIA_BACKENDINFO_HPP

#include <string>

namespace utopia {
    class BackendInfo {
    public:
        const std::string &getName() const {
            return name_;
        }

        void setName(const std::string &name) {
            this->name_ = name;
        }

    private:
        std::string name_;
    };
}

#endif //UTOPIA_UTOPIA_BACKENDINFO_HPP

#ifndef UTOPIA_UTOPIA_BACKENDINFO_HPP
#define UTOPIA_UTOPIA_BACKENDINFO_HPP

#include <string>
#include <utility>
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {
    class BackendInfo {
    public:
        const std::string &get_name() const { return name_; }

        void set_name(const std::string &name) { this->name_ = name; }

        BackendInfo(std::string name = "undefined") : name_(std::move(name)) {}

    private:
        std::string name_;
    };

    template <class T>
    const BackendInfo &backend_info(const T &) {
        return Traits<T>::backend_info();
    }
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_BACKENDINFO_HPP

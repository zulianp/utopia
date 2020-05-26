#ifndef UTOPIA_MAKE_UNIQUE_HPP
#define UTOPIA_MAKE_UNIQUE_HPP

#include <memory>
#include <utility>

namespace utopia {

#ifdef WITH_CPP14

    using std::make_unique;

#else
    // note: this implementation does not disable this overload for array types
    template <typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args &&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

#endif  // WITH_CPP14

    template <typename T>
    inline std::shared_ptr<T> unique_to_shared(std::unique_ptr<T> &ptr) {
        return std::move(ptr);
    }

    template <typename T>
    inline std::shared_ptr<T> unique_to_shared(std::unique_ptr<T> &&ptr) {
        return std::move(ptr);
    }
}  // namespace utopia

#endif  // UTOPIA_MAKE_UNIQUE_HPP

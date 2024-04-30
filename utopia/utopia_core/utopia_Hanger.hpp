#ifndef UTOPIA_HANGER_HPP
#define UTOPIA_HANGER_HPP

#include <map>
#include <memory>
#include <string>

namespace utopia {
    class Hanger {
    public:
        using Key = std::string;
        using Value = std::shared_ptr<void>;

        inline Value &operator[](const Key &k) { return data_[k]; }

        template <typename T>
        inline std::shared_ptr<T> find(const Key &k) const {
            auto it = data_.find(k);
            if (it == data_.end()) return nullptr;
            return std::dynamic_pointer_cast<T>(it->second);
        }

    private:
        std::map<Key, Value> data_;
    };
}  // namespace utopia

#endif  // UTOPIA_HANGER_HPP

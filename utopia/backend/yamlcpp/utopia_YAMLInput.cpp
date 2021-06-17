#include "utopia_YAMLInput.hpp"

#include "utopia_Utils.hpp"
#include "utopia_make_unique.hpp"

#include <yaml-cpp/yaml.h>

namespace utopia {

    class YAMLInput::Impl {
    public:
        YAML::Node root_;
        YAML::Node &root() { return root_; }
        const YAML::Node &root() const { return root_; }
    };

    YAMLInput::YAMLInput(const YAML::Node &node) : impl_(utopia::make_unique<Impl>()) { impl_->root_ = node; }
    YAMLInput::YAMLInput(YAML::Node &&node) : impl_(utopia::make_unique<Impl>()) { impl_->root_ = std::move(node); }

    YAMLInput::YAMLInput() : impl_(utopia::make_unique<Impl>()) {}
    YAMLInput::~YAMLInput() {}

    bool YAMLInput::key_exists(const std::string &key) const {
        auto &node = impl_->root();
        if (node[key]) {
            return true;
        } else {
            return false;
        }
    }

    bool YAMLInput::open(const Path &path) {
        impl_->root_ = YAML::LoadFile(path.to_string());
        return good();
    }

    SizeType YAMLInput::size() const {
        if (good()) {
            return impl_->root().size();
        } else {
            return 0;
        }
    }

    void YAMLInput::get(std::vector<std::shared_ptr<IConvertible>> &values) {
        auto &node = impl_->root();

        for (std::size_t i = 0; i < node.size(); i++) {
            values.push_back(std::make_shared<Convertible<std::string>>(node[i].as<std::string>()));
        }
    }

    void YAMLInput::get_all(std::function<void(Input &)> lambda) {
        auto &node = impl_->root();

        for (std::size_t i = 0; i < node.size(); i++) {
            YAMLInput child(node[i]);
            lambda(child);
        }
    }

    void YAMLInput::get(const std::string &key, bool &val) {
        auto &node = impl_->root();
        if (node[key]) {
            val = node[key].as<bool>();
        }
    }

    void YAMLInput::get(const std::string &key, double &val) {
        auto &node = impl_->root();
        if (node[key]) {
            val = node[key].as<double>();
        }
    }

    void YAMLInput::get(const std::string &key, int &val) {
        auto &node = impl_->root();
        if (node[key]) {
            val = node[key].as<int>();
        }
    }

    void YAMLInput::get(const std::string &key, long &val) {
        auto &node = impl_->root();
        if (node[key]) {
            val = node[key].as<long>();
        }
    }

    void YAMLInput::get(const std::string &key, unsigned long &val) {
        auto &node = impl_->root();
        if (node[key]) {
            val = node[key].as<unsigned long>();
        }
    }

    void YAMLInput::get(const std::string &key, long long int &val) {
        auto &node = impl_->root();
        if (node[key]) {
            val = node[key].as<long long int>();
        }
    }

    void YAMLInput::get(const std::string &key, std::string &val) {
        auto &node = impl_->root();
        if (node[key]) {
            val = node[key].as<std::string>();
        }
    }

    void YAMLInput::get(const std::string &key, Configurable &val) {
        auto &node = impl_->root();
        if (node[key]) {
            YAMLInput branch(node[key]);
            val.read(branch);
        }
    }

    void YAMLInput::get(const std::string &key, std::function<void(Input &)> lambda) {
        auto &node = impl_->root();
        if (node[key]) {
            YAMLInput branch(node[key]);
            lambda(branch);
        }
    }

    bool YAMLInput::good() const { return static_cast<bool>(impl_->root()); }

    bool YAMLInput::is_collection() const {
        auto &node = impl_->root();

        return (node.Type() == YAML::NodeType::Sequence);
    }

}  // namespace utopia
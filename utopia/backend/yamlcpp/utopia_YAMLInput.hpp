#ifndef UTOPIA_YAML_INPUT_HPP
#define UTOPIA_YAML_INPUT_HPP

#include <string>
#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_YAML_CPP

#include <memory>
#include "utopia_Input.hpp"
#include "utopia_Path.hpp"

namespace YAML {
    class Node;
}

namespace utopia {

    class YAMLInput final : public Input {
    public:
        using Input::get;

        YAMLInput(const YAML::Node &node);
        YAMLInput(YAML::Node &&node);

        static std::shared_ptr<YAMLInput> FromString(const std::string &str);

        YAMLInput();
        ~YAMLInput() override;

        bool open(const Path &path);

        SizeType size() const override;

        void get(const std::string &key, std::vector<std::string> &v) override;
        void get(const std::string &key, std::vector<float> &v) override;
        void get(const std::string &key, std::vector<double> &v) override;

        void get(std::vector<std::shared_ptr<IConvertible>> &values) override;
        void get_all(std::function<void(Input &)> lambda) override;
        void get(const std::string &key, bool &val) override;
        void get(const std::string &key, double &val) override;
        void get(const std::string &key, float &val) override;
        void get(const std::string &key, int &val) override;
        void get(const std::string &key, long &val) override;
        void get(const std::string &key, unsigned long &val) override;

        void get(const std::string &key, long long int &val) override;
        // void get(const std::string &key, long long &val) override;

        // void get(const std::string &key, SizeType &val) override;
        void get(const std::string &key, std::string &val) override;
        void get(const std::string &key, Configurable &val) override;
        void get(const std::string &key, std::function<void(Input &)> lambda) override;
        bool good() const override;

        bool key_exists(const std::string &key) const override;

        bool is_collection() const override;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };
}  // namespace utopia

#endif  // UTOPIA_WITH_YAML_CPP
#endif  // UTOPIA_YAML_INPUT_HPP

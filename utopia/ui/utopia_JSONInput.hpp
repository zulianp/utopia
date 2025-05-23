#ifndef UTOPIA_JSON_STREAM_HPP
#define UTOPIA_JSON_STREAM_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_ENABLE_JSON

#include <memory>
#include "utopia_Input.hpp"
#include "utopia_Path.hpp"

namespace utopia {

    class JSONInput final : public Input {
    public:
        using Input::get;

        JSONInput();
        ~JSONInput() override;

        bool open(const Path &path);

        SizeType size() const override;
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

#endif  // UTOPIA_ENABLE_JSON
#endif  // UTOPIA_JSON_STREAM_HPP

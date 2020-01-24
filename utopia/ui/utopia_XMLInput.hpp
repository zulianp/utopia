#ifndef UTOPIA_XML_STREAM_HPP
#define UTOPIA_XML_STREAM_HPP

#include <memory>
#include "utopia_Base.hpp"
#include "utopia_Path.hpp"
#include "utopia_Input.hpp"


namespace utopia {

    class XMLInput final : public Input {
    public:
        XMLInput();
        ~XMLInput();

        bool open(const Path &path);// override;

        void get(bool &val);// override;
        void get(double &val);// override;
        void get(int &val);// override;
        void get(long &val);// override;
        void get(long long &val);// override;
        void get(unsigned long &val);// override;
        // void get(SizeType &val);// override;
        void get(std::string &val);// override;
        void get(Configurable &val);// override;
        void get(std::function<void(Input &)> lambda);// override;

        void get(const std::string &key, bool &val) override;
        void get(const std::string &key, double &val) override;
        void get(const std::string &key, int &val) override;
        void get(const std::string &key, long &val) override;
        void get(const std::string &key, long long &val) override;
        void get(const std::string &key, unsigned long &val) override;
        // void get(const std::string &key, SizeType &val) override;
        void get(const std::string &key, std::string &val) override;
        void get(const std::string &key, Configurable &val) override;
        void get(const std::string &key, std::function<void(Input &)> lambda) override;

        void get_all(std::function<void(Input &)> lambda) override;
        void get(std::vector<std::shared_ptr<IConvertible>> &/*values*/) override {
            assert(false && "implement me");
        }

        bool good() const override;

        SizeType size() const override;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;

        bool object_begin(const std::string &name);
        bool object_end();

        void next();   //override;
        void array_start();  //override;
        void array_finish(); //override;
    };
}

#endif //UTOPIA_XML_STREAM_HPP

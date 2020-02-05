#ifndef UTOPIA_INPUT_STREAM_HPP
#define UTOPIA_INPUT_STREAM_HPP

#include <memory>
#include <string>
#include <vector>
#include <set>
#include <functional>
#include <ostream>

#include "utopia_Base.hpp"
#include "utopia_Convertible.hpp"
#include "utopia_Path.hpp"

namespace utopia {
    class Path;
    class Input;

    enum VerbosityLevel  {  VERBOSITY_LEVEL_QUIET         =-1,
                            VERBOSITY_LEVEL_NORMAL        = 0,
                            VERBOSITY_LEVEL_VERY_VERBOSE  = 1,
                            VERBOSITY_LEVEL_DEBUG         = 2 };

    // enum NormSchedule{  EVERY_ITER = 1,
    //                     NEVER = 2};

    enum MultilevelNormSchedule{    ALWAYS = 1,
                                    OUTER_CYCLE = 2};

    class Configurable {
    public:
        virtual ~Configurable() {}
        virtual void read(Input &is) = 0;
        virtual void print_usage(std::ostream &os = std::cout) const;
        virtual void print_param_usage(std::ostream &os, const std::string & name, const std::string & type, const std::string & description, const std::string & default_settings) const;
        virtual bool import(const Path &path);
        virtual bool import(
            const std::string &key,
            const Path &path);
    };

    class Input  /* : public Clonable */ {
    public:
        Input() {}
        virtual ~Input() {}

        virtual SizeType size() const = 0;
        virtual void get(std::vector<std::shared_ptr<IConvertible>> &values) = 0;
        virtual void get_all(std::function<void(Input &)> lambda) = 0;

        virtual void get(const std::string &key, bool &val) = 0;
        virtual void get(const std::string &key, double &val) = 0;
        virtual void get(const std::string &key, int &val) = 0;
        virtual void get(const std::string &key, long &val) = 0;
        virtual void get(const std::string &key, unsigned long &val) = 0;

        virtual void get(const std::string &key, long long int &val) = 0;
        // virtual void get(const std::string &key, long long &val) = 0;

        // virtual void get(const std::string &key, SizeType &val) = 0;
        virtual void get(const std::string &key, std::string &val) = 0;
        virtual void get(const std::string &key, Configurable &val) = 0;
        virtual void get(const std::string &key, std::function<void(Input &)> lambda) = 0;
        virtual bool good() const = 0;

    private:
        Input(const Input &) {}
        Input &operator=(const Input &) {
            return *this;
        }
    };
}

#endif //UTOPIA_INPUT_STREAM_HPP

#ifndef UTOPIA_INPUT_STREAM_HPP
#define UTOPIA_INPUT_STREAM_HPP

#include <functional>
#include <memory>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "utopia_Base.hpp"
#include "utopia_Convertible.hpp"
#include "utopia_Path.hpp"

namespace utopia {
    class Path;
    class Input;

    enum VerbosityLevel {
        VERBOSITY_LEVEL_QUIET = -1,
        VERBOSITY_LEVEL_NORMAL = 0,
        VERBOSITY_LEVEL_VERY_VERBOSE = 1,
        VERBOSITY_LEVEL_DEBUG = 2
    };

    // enum NormSchedule{  EVERY_ITER = 1,
    //                     NEVER = 2};

    enum MultilevelNormSchedule { ALWAYS = 1, OUTER_CYCLE = 2 };

    class Configurable {
    public:
        virtual ~Configurable() = default;
        virtual void read(Input &is) = 0;
        virtual void print_usage(std::ostream &os = std::cout) const;
        virtual void print_param_usage(std::ostream &os,
                                       const std::string &name,
                                       const std::string &type,
                                       const std::string &description,
                                       const std::string &default_settings) const;
        virtual bool import(const Path &path);
        virtual bool import(const std::string &key, const Path &path);
    };

    class Input /* : public Clonable */ {
    public:
        Input() = default;
        virtual ~Input() = default;

        virtual SizeType size() const = 0;
        virtual void get(std::vector<std::shared_ptr<IConvertible>> &values) = 0;
        virtual void get_all(std::function<void(Input &)> lambda) = 0;
        virtual bool is_collection() const = 0;

        virtual void get(const std::string &key, Path &val) { get(key, val.raw_type()); }
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

        virtual bool key_exists(const std::string &key) const = 0;

        class Attribute {
        public:
            virtual ~Attribute() = default;

            virtual void key_exists() const {}
            virtual void key_does_not_exist() const {}
        };

        class KeyDeprecated final : public Attribute {
        public:
            KeyDeprecated(std::string deprecated_key, std::string new_key)
                : deprecated_key_(std::move(deprecated_key)), new_key_(std::move(new_key)) {}

            void key_exists() const override;
            std::string deprecated_key_;
            std::string new_key_;
        };

        class KeyRequired final : public Attribute {
        public:
            KeyRequired(std::string key) : key_(std::move(key)) {}

            void key_does_not_exist() const override;

            std::string key_;
        };

        template <typename T>
        void get(const std::string &key, T &val, const Attribute &attribute) {
            if (key_exists(key)) {
                attribute.key_exists();
                get(key, val);
            } else {
                attribute.key_does_not_exist();
            }
        }

        template <typename T>
        void get_deprecated(const std::string &old_key, const std::string &new_key, T &val) {
            get(old_key, val, KeyDeprecated(old_key, new_key));
        }

        template <typename T>
        void require(const std::string &key, T &val) {
            get(key, val, KeyRequired(key));
        }

        void require(const std::string &key, std::function<void(Input &)> lambda) {
            get(key, lambda, KeyRequired(key));
        }

        virtual bool good() const = 0;

        Input &operator=(Input &&) = default;
        Input(Input &&) = default;

    private:
        Input(const Input &) = default;
        Input &operator=(const Input &) = default;
    };
}  // namespace utopia

#endif  // UTOPIA_INPUT_STREAM_HPP

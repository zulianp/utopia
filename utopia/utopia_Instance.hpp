#ifndef UTOPIA_UTOPIA_INSTANCE_HPP
#define UTOPIA_UTOPIA_INSTANCE_HPP

#include "utopia_Input.hpp"
#include "utopia_Library.hpp"
#include "utopia_Logger.hpp"

#include <cassert>
#include <list>
#include <map>
#include <string>
#include <vector>

namespace utopia {
    class Utopia final : public Configurable {
    public:
        static void Init(int argc, char *argv[]);
        static int Finalize();
        static void Abort();
        static void Abort(const std::string &message);

        void read(Input &is) override;
        void print_usage(std::ostream &os) const override;

        inline std::string get(const std::string &key) const {
            auto it = settings_.find(key);
            if (it == settings_.end()) return "";
            return it->second;
        }

        inline void set(const std::string &key, const std::string &value) { settings_[key] = value; }

        static Utopia &instance();

        bool verbose() const;

        Logger &logger() {
            assert(logger_);
            return *logger_;
        }

        Logger &maintenance_logger() {
            assert(logger_);
            return *maintenance_logger_;
        }

        void set_exit_code(const int code) { exit_code_ = code; }

        void read_input(int argc, char *argv[]);

        inline void add_library(std::unique_ptr<Library> &&l) { libraries_.push_back(std::move(l)); }

    private:
        Utopia();
        std::map<std::string, std::string> settings_;
        std::shared_ptr<Logger> logger_;
        std::shared_ptr<Logger> maintenance_logger_;
        std::list<std::unique_ptr<Library>> libraries_;
        int exit_code_{EXIT_SUCCESS};

        inline void add_library_with_priority(std::unique_ptr<Library> &&l) { libraries_.push_front(std::move(l)); }
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_INSTANCE_HPP

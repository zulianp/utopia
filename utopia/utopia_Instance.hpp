#ifndef UTOPIA_UTOPIA_INSTANCE_HPP
#define UTOPIA_UTOPIA_INSTANCE_HPP

#include "utopia_Logger.hpp"

#include <map>
#include <string>
#include <cassert>

namespace utopia {
    class Utopia final {
    public:
        static void Init(int argc, char *argv[]);
        static int Finalize();
        
        inline std::string get(const std::string &key) const
        {
            auto it = settings_.find(key);
            if(it == settings_.end()) return "";
            return it->second;
        }
        
        inline void set(const std::string &key, const std::string &value)
        {
            settings_[key] = value;
        }
        
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

        void set_exit_code(const int code)
        {
            exit_code_ = code;
        }

    private:
        Utopia();
        std::map<std::string, std::string> settings_;
        std::shared_ptr<Logger> logger_;
        std::shared_ptr<Logger> maintenance_logger_;
        int exit_code_;
    };
}

#endif //UTOPIA_UTOPIA_INSTANCE_HPP

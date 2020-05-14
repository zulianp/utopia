#ifndef UTOPIA_LOGGER_HPP
#define UTOPIA_LOGGER_HPP

#include <string>
#include <iostream>
#include <memory>
#include <vector>

#define utopia_warning(macro_msg_) { utopia::Utopia::instance().logger().warning(__FILE__, __LINE__, macro_msg_); }
#define utopia_error(macro_msg_)   { utopia::Utopia::instance().logger().error(__FILE__, __LINE__, macro_msg_);   }
#define utopia_status(macro_msg_)  { utopia::Utopia::instance().logger().status(__FILE__, __LINE__, macro_msg_);  }

//maintener messages (the once version is printed only one time per run)
#define m_utopia_warning(macro_msg_) { utopia::Utopia::instance().maintenance_logger().warning(__FILE__, __LINE__, macro_msg_); }
#define m_utopia_warning_once(macro_msg_) { static bool msg_printed_ = false; if(!msg_printed_) { m_utopia_warning(macro_msg_); msg_printed_ = true; }}

#define m_utopia_error(macro_msg_)   { utopia::Utopia::instance().maintenance_logger().error(__FILE__, __LINE__, macro_msg_);   }
#define m_utopia_error_once(macro_msg_) { static bool msg_printed_ = false; if(!msg_printed_) { m_utopia_error(macro_msg_); msg_printed_ = true; }}

#define m_utopia_status(macro_msg_)  { utopia::Utopia::instance().maintenance_logger().status(__FILE__, __LINE__, macro_msg_);  }
#define m_utopia_status_once(macro_msg_) { static bool msg_printed_ = false; if(!msg_printed_) { m_utopia_status(macro_msg_); msg_printed_ = true; }}

namespace utopia {

    class Logger {
    public:
        class Entry;

        virtual void warning(
            const std::string &file,
            const int line,
            const std::string &message) = 0;

        virtual void error(
            const std::string &file,
            const int line,
            const std::string &message) = 0;

        virtual void status(
            const std::string &file,
            const int line,
            const std::string &message) = 0;

        virtual void flush() = 0;

        Logger();
        virtual ~Logger();
    };

    class NullLogger final : public Logger {
    public:
        inline void warning(
            const std::string &,
            const int,
            const std::string &) override {}

        inline void error(
            const std::string &,
            const int,
            const std::string &) override {}

        inline void status(
            const std::string &,
            const int,
            const std::string &) override {}

        inline void flush() override {}
    };

    class StandardLogger final : public Logger {
    public:
        void status(
            const std::string &file,
            const int line,
            const std::string &message) override;

        void warning(
            const std::string &file,
            const int line,
            const std::string &message) override;

        void error(
            const std::string &file,
            const int line,
            const std::string &message) override;

        StandardLogger(
            std::ostream &status_stream = std::cout,
            std::ostream &warning_stream = std::cerr,
            std::ostream &error_stream = std::cerr
            );

        ~StandardLogger() override;

        void flush() override;

        inline void set_direct_output(
            const bool status,
            const bool warning,
            const bool error)
        {
            status_direct_output_ = status;
            warning_direct_output_ = warning;
            error_direct_output_ = error;
        }

    private:
        std::ostream &status_stream_;
        std::ostream &warning_stream_;
        std::ostream &error_stream_;

        bool status_direct_output_{true};
        bool warning_direct_output_{true};
        bool error_direct_output_{true};

        std::vector<std::shared_ptr<Logger::Entry>> status_entries_;
        std::vector<std::shared_ptr<Logger::Entry>> warning_entries_;
        std::vector<std::shared_ptr<Logger::Entry>> error_entries_;
    };
}

#endif //UTOPIA_LOGGER_HPP

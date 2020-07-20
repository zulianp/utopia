#include "utopia_Logger.hpp"

#include <utility>

namespace utopia {
    Logger::Logger() = default;
    Logger::~Logger() = default;

    class Logger::Entry final {
    public:
        inline explicit Entry(std::string file, const int line, std::string message)
            : file(std::move(file)), line(line), message(std::move(message)) {}

        inline friend std::ostream &operator<<(std::ostream &os, const Entry &e) {
            os << e.file << ":" << e.line << "\n" << e.message;
            return os;
        }

        std::string file;
        int line;
        std::string message;
    };

    void StandardLogger::status(const std::string &file, const int line, const std::string &message) {
        if (status_direct_output_) {
            status_stream_ << "[Status] " << Entry(file, line, message) << std::endl;
        } else {
            status_entries_.push_back(std::make_shared<Entry>(file, line, message));
        }
    }

    void StandardLogger::warning(const std::string &file, const int line, const std::string &message) {
        if (warning_direct_output_) {
            warning_stream_ << "[Warning] " << Entry(file, line, message) << std::endl;
        } else {
            warning_entries_.push_back(std::make_shared<Entry>(file, line, message));
        }
    }

    void StandardLogger::error(const std::string &file, const int line, const std::string &message) {
        if (error_direct_output_) {
            error_stream_ << "[Error] " << Entry(file, line, message) << std::endl;
        } else {
            error_entries_.push_back(std::make_shared<Entry>(file, line, message));
        }
    }

    void StandardLogger::flush() {
        if (!status_entries_.empty()) {
            status_stream_ << "--------------------------------------------------------" << std::endl;
            status_stream_ << "[Status]: " << std::endl;
            for (auto e : status_entries_) {
                status_stream_ << *e << std::endl;
            }
            status_stream_ << "--------------------------------------------------------" << std::endl;
            status_stream_.clear();
        }

        if (!warning_entries_.empty()) {
            status_stream_ << "--------------------------------------------------------" << std::endl;
            warning_stream_ << "[Warning]:\n" << std::endl;
            for (auto e : warning_entries_) {
                warning_stream_ << *e << "\n" << std::endl;
            }

            status_stream_ << "--------------------------------------------------------" << std::endl;
            warning_stream_.clear();
        }
        if (!error_entries_.empty()) {
            status_stream_ << "--------------------------------------------------------" << std::endl;
            error_stream_ << "[Error]: " << std::endl;
            for (auto e : error_entries_) {
                error_stream_ << *e << std::endl;
            }

            status_stream_ << "--------------------------------------------------------" << std::endl;
            error_entries_.clear();
        }
    }

    StandardLogger::StandardLogger(std::ostream &status_stream,
                                   std::ostream &warning_stream,
                                   std::ostream &error_stream)
        : status_stream_(status_stream),
          warning_stream_(warning_stream),
          error_stream_(error_stream)

    {}

    StandardLogger::~StandardLogger() = default;
}  // namespace utopia

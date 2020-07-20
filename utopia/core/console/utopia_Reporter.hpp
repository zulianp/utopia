#ifndef UTOPIA_REPORTER_HPP
#define UTOPIA_REPORTER_HPP

#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>

#include "utopia_Input.hpp"
#include "utopia_Utils.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_IOStream.hpp"

namespace utopia {

    class NullOStream final : public OStream {
    public:
        inline void write(bool) override {}
        inline void write(short) override {}
        inline void write(unsigned short) override {}
        inline void write(int) override {}
        inline void write(unsigned int) override {}
        inline void write(long) override {}
        inline void write(long long) override {}
        inline void write(unsigned long) override {}
        inline void write(unsigned long long) override{};
        inline void write(float) override {}
        inline void write(double) override {}
        inline void write(long double) override {}
        inline void write(void *) override {}

        // stream buffers (2)
        inline void write(std::streambuf *) override {}
        // manipulators (3)
        inline void write(std::ostream &(*)(std::ostream &)) override {}
        inline void write(std::ios &(*)(std::ios &)) override {}
        inline void write(std::ios_base &(*)(std::ios_base &)) override {}

        // single character (1)
        inline void write(char) override {}
        inline void write(signed char) override {}
        inline void write(unsigned char) override {}
        // character sequence (2)
        inline void write(const char *) override {}
        inline void write(const signed char *) override {}
        inline void write(const unsigned char *) override {}

        inline void write(const Describable &) override {}
    };

    class StreamWrapper : public OStream {
    public:
        ~StreamWrapper() = default;

        inline void write(bool val) override { stream() << val; }
        inline void write(short val) override { stream() << val; }
        inline void write(unsigned short val) override { stream() << val; }
        inline void write(int val) override { stream() << val; }
        inline void write(unsigned int val) override { stream() << val; }
        inline void write(long val) override { stream() << val; }
        inline void write(long long val) override { stream() << val; }
        inline void write(unsigned long val) override { stream() << val; }
        inline void write(unsigned long long val) override { stream() << val; };
        inline void write(float val) override { stream() << val; }
        inline void write(double val) override { stream() << val; }
        inline void write(long double val) override { stream() << val; }
        inline void write(void *val) override { stream() << val; }

        // stream buffers (2)
        inline void write(std::streambuf *sb) override { stream() << sb; }
        // manipulators (3)
        inline void write(std::ostream &(*pf)(std::ostream &)) override { stream() << pf; }
        inline void write(std::ios &(*pf)(std::ios &)) override { stream() << pf; }
        inline void write(std::ios_base &(*pf)(std::ios_base &)) override { stream() << pf; }

        // single character (1)
        inline void write(char c) override { stream() << c; }
        inline void write(signed char c) override { stream() << c; }
        inline void write(unsigned char c) override { stream() << c; }
        // character sequence (2)
        inline void write(const char *s) override { stream() << s; }
        inline void write(const signed char *s) override { stream() << s; }
        inline void write(const unsigned char *s) override { stream() << s; }

        inline void write(const Describable &d) override { d.describe(stream()); }

        inline std::ostream &stream() { return *stream_ptr_; }

        inline void wrap_cout() { stream_ptr_ = &std::cout; }
        inline void wrap_cerr() { stream_ptr_ = &std::cerr; }
        inline void wrap(std::ostream &os) { stream_ptr_ = &os; }

        StreamWrapper() : stream_ptr_(&std::cout) {}
        StreamWrapper(std::ostream &os) : stream_ptr_(&os) {}

    private:
        std::ostream *stream_ptr_;
    };

    class FileLogger final : public StreamWrapper {
    public:
        ~FileLogger() {
            os_ << std::flush;
            os_.close();
        }

        bool open(const Path &path) {
            os_.open(path.c_str());
            this->wrap(os_);
            return !os_.good();
        }

    private:
        std::ofstream os_;
    };

    class Reporter : public Configurable {
    public:
        inline static Reporter &instance() {
            static Reporter instance_;
            return instance_;
        }

        ~Reporter() = default;

        void read(Input &in) override;

        inline OStream &cout() { return *cout_; }
        inline OStream &cerr() { return *cerr_; }
        inline OStream &dev() { return *dev_; }

        Reporter()
            : cout_(utopia::make_unique<StreamWrapper>(std::cout)),
              cerr_(utopia::make_unique<StreamWrapper>(std::cerr)),
              dev_(utopia::make_unique<StreamWrapper>(std::cout)) {}

    private:
        std::unique_ptr<OStream> cout_;
        std::unique_ptr<OStream> cerr_;
        std::unique_ptr<OStream> dev_;
    };

}  // namespace utopia

#endif  // UTOPIA_REPORTER_HPP

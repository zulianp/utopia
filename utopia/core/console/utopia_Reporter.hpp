#ifndef UTOPIA_REPORTER_HPP
#define UTOPIA_REPORTER_HPP

#include <fstream>
#include <iostream>
#include <ostream>

#include "utopia_Input.hpp"
#include "utopia_Utils.hpp"
#include "utopia_make_unique.hpp"

namespace utopia {

    class AppOutputStream : public Configurable {
    public:
        virtual ~AppOutputStream() = default;
        AppOutputStream() = default;

        AppOutputStream(const AppOutputStream &) = delete;
        AppOutputStream &operator=(const AppOutputStream &) = delete;

        template <typename T>
        inline AppOutputStream &operator<<(const T &val) {
            this->write(val);
            return *this;
        }

        virtual void write(bool val) = 0;
        virtual void write(short val) = 0;
        virtual void write(unsigned short val) = 0;
        virtual void write(int val) = 0;
        virtual void write(unsigned int val) = 0;
        virtual void write(long val) = 0;
        virtual void write(unsigned long val) = 0;
        virtual void write(float val) = 0;
        virtual void write(double val) = 0;
        virtual void write(long double val) = 0;
        virtual void write(void *val) = 0;

        // stream buffers (2)
        virtual void write(std::streambuf *sb) = 0;
        // manipulators (3)
        virtual void write(std::ostream &(*pf)(std::ostream &)) = 0;
        virtual void write(std::ios &(*pf)(std::ios &)) = 0;
        virtual void write(std::ios_base &(*pf)(std::ios_base &)) = 0;

        // single character (1)
        virtual void write(char c) = 0;
        virtual void write(signed char c) = 0;
        virtual void write(unsigned char c) = 0;
        // character sequence (2)
        virtual void write(const char *s) = 0;
        virtual void write(const signed char *s) = 0;
        virtual void write(const unsigned char *s) = 0;

        virtual void read(Input &) {}
    };

    class NullAppOutputStream final : public AppOutputStream {
    public:
        inline void write(bool) override {}
        inline void write(short) override {}
        inline void write(unsigned short) override {}
        inline void write(int) override {}
        inline void write(unsigned int) override {}
        inline void write(long) override {}
        inline void write(unsigned long) override {}
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
    };

    class StreamWrapper final : public AppOutputStream {
    public:
        template <typename T>
        struct CloseFile {
            CloseFile() /* noexcept */
                = default;

            template <typename U>
            CloseFile(const CloseFile<U> &,
                      typename std::enable_if<std::is_convertible<U *, T *>::value>::type * = nullptr) noexcept {}

            void operator()(T *const obj) const /* noexcept */
            {
                // do nothing
                obj->close();
                delete obj;
            }
        };

        ~StreamWrapper() = default;

        inline void write(bool val) override { stream() << val; }
        inline void write(short val) override { stream() << val; }
        inline void write(unsigned short val) override { stream() << val; }
        inline void write(int val) override { stream() << val; }
        inline void write(unsigned int val) override { stream() << val; }
        inline void write(long val) override { stream() << val; }
        inline void write(unsigned long val) override { stream() << val; }
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

        std::ostream &stream() { return *stream_ptr_; }

        void wrap_cout() { stream_ptr_ = &std::cout; }
        void wrap_cerr() { stream_ptr_ = &std::cerr; }

        StreamWrapper() : stream_ptr_(&std::cout) {}
        StreamWrapper(std::ostream &os) : stream_ptr_(&os) {}

    private:
        std::ostream *stream_ptr_;
    };

    class Reporter : public Configurable {
    public:
        inline static Reporter &instance() {
            static Reporter instance_;
            return instance_;
        }

        ~Reporter() = default;

        void read(Input &) override {}

        inline AppOutputStream &cout() { return *cout_; }
        inline AppOutputStream &cerr() { return *cerr_; }
        inline AppOutputStream &dev() { return *dev_; }

        Reporter()
            : cout_(utopia::make_unique<StreamWrapper>(std::cout)),
              cerr_(utopia::make_unique<StreamWrapper>(std::cerr)),
              dev_(utopia::make_unique<StreamWrapper>(std::cout)) {}

    private:
        std::unique_ptr<AppOutputStream> cout_;
        std::unique_ptr<AppOutputStream> cerr_;
        std::unique_ptr<AppOutputStream> dev_;
    };

    inline AppOutputStream &out() { return Reporter::instance().cout(); }
    inline AppOutputStream &err() { return Reporter::instance().cerr(); }
    inline AppOutputStream &dev() { return Reporter::instance().dev(); }

}  // namespace utopia

#endif  // UTOPIA_REPORTER_HPP

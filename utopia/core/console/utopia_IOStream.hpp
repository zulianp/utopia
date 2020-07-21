#ifndef UTOPIA_IO_STREAM_HPP
#define UTOPIA_IO_STREAM_HPP

#include <ostream>
#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"

namespace utopia {

    class OStream : public Configurable {
    public:
        virtual ~OStream() = default;
        OStream() = default;

        OStream(const OStream &) = delete;
        OStream &operator=(const OStream &) = delete;

        inline OStream &operator<<(std::ostream &(*pf)(std::ostream &)) {
            write(pf);
            return *this;
        }

        inline OStream &operator<<(std::ios &(*pf)(std::ios &)) {
            write(pf);
            return *this;
        }

        inline OStream &operator<<(std::ios_base &(*pf)(std::ios_base &)) {
            write(pf);
            return *this;
        }

        inline OStream &operator<<(const Describable &val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(bool val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(short val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(unsigned short val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(int val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(unsigned int val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(long val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(long long val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(unsigned long val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(unsigned long long val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(float val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(double val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(long double val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(void *val) {
            write(val);
            return *this;
        }

        inline OStream &operator<<(const std::string &val) {
            write(val);
            return *this;
        }

        // template <typename T>
        // inline OStream &operator<<(const T &val) {
        //     this->write(val);
        //     return *this;
        // }

        // inline OStream &operator<<(std::ostream &(*pf)(std::ostream &)) {
        //     write(pf);
        //     return *this;
        // }

        // inline OStream &operator<<(std::ios &(*pf)(std::ios &)) {
        //     write(pf);
        //     return *this;
        // }

        // inline OStream &operator<<(std::ios_base &(*pf)(std::ios_base &)) {
        //     write(pf);
        //     return *this;
        // }

        virtual void write(const Describable &val) = 0;
        virtual void write(bool val) = 0;
        virtual void write(short val) = 0;
        virtual void write(unsigned short val) = 0;
        virtual void write(int val) = 0;
        virtual void write(unsigned int val) = 0;
        virtual void write(long val) = 0;
        virtual void write(long long val) = 0;
        virtual void write(unsigned long val) = 0;
        virtual void write(unsigned long long val) = 0;
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
        virtual void write(const std::string &s) { write(s.c_str()); }
        virtual void write(const signed char *s) = 0;
        virtual void write(const unsigned char *s) = 0;

        virtual void read(Input &) {}
    };

    // we only use these directly
    OStream &out();
    OStream &err();
    OStream &dev();
}  // namespace utopia

#endif  // UTOPIA_IO_STREAM_HPP

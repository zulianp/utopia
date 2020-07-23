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

        template <typename T>
        OStream &operator<<(const T &v) {
            if (!muted_) {
                stream() << v;
            }

            return *this;
        }

        inline OStream &operator<<(std::ostream &(*pf)(std::ostream &)) {
            if (!muted_) {
                stream() << pf;
            }
            return *this;
        }

        inline OStream &operator<<(std::ios &(*pf)(std::ios &)) {
            if (!muted_) {
                stream() << pf;
            }
            return *this;
        }

        inline OStream &operator<<(std::ios_base &(*pf)(std::ios_base &)) {
            if (!muted_) {
                stream() << pf;
            }
            return *this;
        }

        virtual std::ostream &stream() = 0;
        virtual void read(Input &) {}

    protected:
        bool muted_{false};
    };

    // we only use these directly
    OStream &out();
    OStream &err();
    OStream &dev();
}  // namespace utopia

#endif  // UTOPIA_IO_STREAM_HPP

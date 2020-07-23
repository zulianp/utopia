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

    class StreamWrapper : public OStream {
    public:
        ~StreamWrapper() = default;

        inline std::ostream &stream() override { return *stream_ptr_; }

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

    class NullOStream final : public StreamWrapper {
    public:
        NullOStream() { this->muted_ = true; }
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

#ifndef UTOPIA_CL_TRAITS_HPP
#define UTOPIA_CL_TRAITS_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

#include "cl.hpp"
#include "utopia_opencl_Error.hpp"

namespace utopia {

    class CLContext {
    public:
        inline cl::Device &current_device() { return devices_[current_device_]; }

        inline cl::Platform &current_platform() { return platforms_[current_platform_]; }

        inline const cl::Device &current_device() const { return devices_[current_device_]; }

        inline const cl::Platform &current_platform() const { return platforms_[current_platform_]; }

        inline cl::CommandQueue &current_queue() { return queue_; }

        inline cl::Context &current() { return current_; }

        inline static CLContext &instance() {
            static CLContext instance_;
            return instance_;
        }

        void describe_all_resources() const {
            utopia::out() << "platforms:\n";

            for (const auto &p : platforms_) {
                utopia::out() << "----------------------------------------------\n";
                utopia::out() << "[" << p.getInfo<CL_PLATFORM_NAME>() << "] devices: "
                              << "\n";

                std::vector<cl::Device> all_devices;
                p.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
                for (const auto &d : all_devices) {
                    utopia::out() << "\t\t- " << d.getInfo<CL_DEVICE_VENDOR>() << " " << d.getInfo<CL_DEVICE_NAME>()
                                  << "\n";
                }

                utopia::out() << "----------------------------------------------\n";
            }
        }

        void describe_current_setup() const {
            utopia::out() << "Using platform: " << current_platform().getInfo<CL_PLATFORM_NAME>() << "\n";
            utopia::out() << "  Using device: " << current_device().getInfo<CL_DEVICE_VENDOR>() << " ";
            utopia::out() << current_device().getInfo<CL_DEVICE_NAME>() << " ";
            utopia::out() << "(" << current_device().getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << " compute units)\n";
        }

        void set_current_device(const long current_device) {
            current_device_ = current_device;
            init();
        }

    private:
        void init() {
            cl::Platform::get(&platforms_);
            current_platform().getDevices(CL_DEVICE_TYPE_ALL, &devices_);
            current_ = cl::Context({current_device()});
            queue_ = cl::CommandQueue(current(), current_device());

            static const bool verbose = false;

            if (verbose) {
                describe_all_resources();
                describe_current_setup();
            }
        }

        CLContext() : current_platform_(0), current_device_(0) { init(); }

        std::vector<cl::Platform> platforms_;
        std::vector<cl::Device> devices_;
        cl::Context current_;

        cl::CommandQueue queue_;

        long current_platform_;
        long current_device_;
    };

    inline bool compile_opencl_program(const std::string &code, cl::Program &program) {
        cl::Program::Sources sources;
        sources.push_back({code.c_str(), code.length()});

        program = cl::Program(CLContext::instance().current(), sources);
        if (program.build({CLContext::instance().current_device()}) != CL_SUCCESS) {
            std::cerr << " Error building: "
                      << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(CLContext::instance().current_device()) << "\n";
            assert(false);
            return false;
        }

        return true;
    }

    inline bool compile_opencl_programs(const std::vector<std::string> &codes, cl::Program &program) {
        cl::Program::Sources sources;

        for (const auto &code : codes) {
            sources.push_back({code.c_str(), code.length()});
        }

        program = cl::Program(CLContext::instance().current(), sources);
        if (program.build({CLContext::instance().current_device()}) != CL_SUCCESS) {
            std::cerr << " Error building: "
                      << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(CLContext::instance().current_device()) << "\n";
            assert(false);
            return false;
        }

        return true;
    }

    inline bool compile_opencl_programs(const cl::Program::Sources &sources,
                                        cl::Program &program,
                                        const std::string &flags) {
        program = cl::Program(CLContext::instance().current(), sources);
        if (program.build({CLContext::instance().current_device()}, flags.c_str()) != CL_SUCCESS) {
            std::cerr << " Error building: "
                      << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(CLContext::instance().current_device()) << "\n";
            assert(false);
            return false;
        }

        return true;
    }

    template <typename T>
    class CLMatrix {
    public:
        SizeType rows, cols;

        std::vector<T> entries;
        cl::Buffer buffer;

        bool resize(const SizeType rows, const SizeType cols) {
            const SizeType n_entries = rows * cols;

            if (n_entries != entries.size()) {
                // entries.reserve(n_entries);
                entries.resize(n_entries);

                this->rows = rows;
                this->cols = cols;

                cl_int err;
                buffer = cl::Buffer(
                    CLContext::instance().current(), CL_MEM_USE_HOST_PTR, sizeof(T) * n_entries, &entries[0], &err);

                assert(check_cl_error(err));
                return CL_SUCCESS == err;
            }

            return true;
        }

        inline bool synch_read_buffer(cl::CommandQueue &queue) {
            cl_int err = queue.enqueueReadBuffer(buffer, CL_TRUE, 0, sizeof(T) * entries.size(), &entries[0]);
            assert(check_cl_error(err));
            return CL_SUCCESS == err;
        }

        inline bool asynch_read_buffer(cl::CommandQueue &queue, cl::Event &event) {
            cl_int err =
                queue.enqueueReadBuffer(buffer, CL_FALSE, 0, sizeof(T) * entries.size(), &entries[0], nullptr, &event);
            assert(check_cl_error(err));
            return CL_SUCCESS == err;
        }

        inline bool synch_write_buffer(cl::CommandQueue &queue) {
            cl_int err = queue.enqueueWriteBuffer(buffer, CL_TRUE, 0, sizeof(T) * entries.size(), &entries[0]);
            assert(check_cl_error(err));
            return CL_SUCCESS == err;
        }

        inline bool asynch_write_buffer(cl::CommandQueue &queue, cl::Event &event) {
            cl_int err =
                queue.enqueueWriteBuffer(buffer, CL_FALSE, 0, sizeof(T) * entries.size(), &entries[0], nullptr, &event);
            assert(check_cl_error(err));
            return CL_SUCCESS == err;
        }
    };

    template <typename T>
    class CLVector {
    public:
        // SizeType size;

        std::vector<T> entries;
        cl::Buffer buffer;

        CLVector() {}

        CLVector(const SizeType n_entries) { resize(n_entries); }

        bool resize(const SizeType n_entries) {
            if (n_entries != entries.size()) {
                entries.reserve(n_entries);
                entries.resize(n_entries);

                cl_int err;
                buffer = cl::Buffer(
                    CLContext::instance().current(), CL_MEM_USE_HOST_PTR, sizeof(T) * n_entries, &entries[0], &err);

                assert(check_cl_error(err));
                return CL_SUCCESS == err;
            }

            return true;
        }

        inline bool synch_read_buffer(cl::CommandQueue &queue) {
            cl_int err = queue.enqueueReadBuffer(buffer, CL_TRUE, 0, sizeof(T) * entries.size(), &entries[0]);
            assert(check_cl_error(err));
            return CL_SUCCESS == err;
        }

        inline bool asynch_read_buffer(cl::CommandQueue &queue, cl::Event &event) {
            cl_int err =
                queue.enqueueReadBuffer(buffer, CL_FALSE, 0, sizeof(T) * entries.size(), &entries[0], nullptr, &event);
            assert(check_cl_error(err));
            return CL_SUCCESS == err;
        }

        inline bool synch_write_buffer(cl::CommandQueue &queue) {
            cl_int err = queue.enqueueWriteBuffer(buffer, CL_TRUE, 0, sizeof(T) * entries.size(), &entries[0]);
            assert(check_cl_error(err));
            return CL_SUCCESS == err;
        }

        inline bool asynch_write_buffer(cl::CommandQueue &queue, cl::Event &event) {
            cl_int err =
                queue.enqueueWriteBuffer(buffer, CL_FALSE, 0, sizeof(T) * entries.size(), &entries[0], nullptr, &event);
            assert(check_cl_error(err));
            return CL_SUCCESS == err;
        }
    };

    template <typename T>
    class CLTraits {
    public:
        typedef T Scalar;
        typedef utopia::CLMatrix<T> Matrix;
        typedef utopia::CLVector<T> Vector;
        typedef utopia::SizeType SizeType;

        enum { Backend = OPENCL_TAG };
    };

    UTOPIA_MAKE_TRAITS_TPL_1(CLVector, CLTraits);
    UTOPIA_MAKE_TRAITS_DENSE_TPL_1(CLMatrix, CLTraits);
    UTOPIA_MAKE_TRAITS(std::shared_ptr<CLVector<double> >, CLTraits<double>);
    UTOPIA_MAKE_TRAITS_DENSE(std::shared_ptr<CLMatrix<double> >, CLTraits<double>);
    UTOPIA_MAKE_TRAITS(std::shared_ptr<CLVector<float> >, CLTraits<float>);
    UTOPIA_MAKE_TRAITS_DENSE(std::shared_ptr<CLMatrix<float> >, CLTraits<float>);

    typedef utopia::Wrapper<utopia::CLMatrix<double>, 2> CLMatrixd;
    typedef utopia::Wrapper<utopia::CLVector<double>, 1> CLVectord;

    template <typename T>
    inline const CLVector<T> &raw_type(const Wrapper<CLVector<T>, 1> &v) {
        return v.implementation();
    }

    template <typename T>
    inline CLVector<T> &raw_type(Wrapper<CLVector<T>, 1> &v) {
        return v.implementation();
    }

    template <typename T>
    inline const CLMatrix<T> &raw_type(const Wrapper<CLMatrix<T>, 2> &v) {
        return v.implementation();
    }

    template <typename T>
    inline CLMatrix<T> &raw_type(Wrapper<CLMatrix<T>, 2> &v) {
        return v.implementation();
    }

    template <typename T>
    inline CLVector<T> &raw_type(const Wrapper<std::shared_ptr<CLVector<T> >, 1> &v) {
        return *v.implementation();
    }

    template <typename T>
    inline CLMatrix<T> &raw_type(const Wrapper<std::shared_ptr<CLMatrix<T> >, 2> &v) {
        return *v.implementation();
    }

    template <typename T>
    inline void disp(Wrapper<CLVector<T>, 1> &w) {
        w.implementation().synch_read_buffer(CLContext::instance().current_queue());
        for (auto v : w.implementation().entries) {
            utopia::out() << v << " ";
        }

        utopia::out() << std::endl;
    }

    template <typename T>
    inline void disp(Wrapper<CLMatrix<T>, 2> &w) {
        w.implementation().synch_read_buffer(CLContext::instance().current_queue());
        auto const &m = w.implementation();
        for (SizeType i = 0; i < m.rows; ++i) {
            const SizeType offset = i * m.cols;
            for (SizeType j = 0; j < m.cols; ++j) {
                utopia::out() << w.implementation().entries[offset + j] << " ";
            }

            utopia::out() << "\n";
        }

        utopia::out() << std::endl;
    }

    class CLStats {
    public:
        inline void describe(std::ostream &os) const {
            os << "-------------------------------------------------\n";
            os << "-------------------- OpenCL stats ---------------\n";
            os << "-------------------------------------------------\n";

            os << "n_kernel_calls:\t\t" << n_kernel_calls << "\n";
            os << "kernel_execution_time:\t" << kernel_execution_time << " seconds\n";

            os << "-------------------------------------------------\n";

            os << "n_code_gen_calls:\t" << n_code_gen_calls << "\n";
            os << "code_generation_time:\t" << code_generation_time << " seconds\n";

            os << "-------------------------------------------------\n";
        }

        void kernel_execution_begin() {
            ++n_kernel_calls;
            chrono_.start();
        }

        void kernel_execution_end() {
            chrono_.stop();
            kernel_execution_time += chrono_.get_seconds();
        }

        void code_generation_begin() {
            ++n_code_gen_calls;
            chrono_.start();
        }

        void code_generation_end() {
            chrono_.stop();
            code_generation_time += chrono_.get_seconds();
        }

        void clear() {
            n_kernel_calls = 0;
            n_code_gen_calls = 0;
            code_generation_time = 0;
            kernel_execution_time = 0;
        }

        static CLStats &instance() {
            static CLStats instance_;
            return instance_;
        }

    private:
        Chrono chrono_;
        long n_kernel_calls;
        long n_code_gen_calls;
        double code_generation_time;
        double kernel_execution_time;

        inline CLStats() : n_kernel_calls(0), n_code_gen_calls(0), code_generation_time(0), kernel_execution_time(0) {}
    };

}  // namespace utopia

#endif  // UTOPIA_CL_TRAITS_HPP

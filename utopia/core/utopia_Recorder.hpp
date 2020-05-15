#ifndef UTOPIA_RECORDER_HPP
#define UTOPIA_RECORDER_HPP

#include <iomanip>
#include <ostream>
#include "utopia_Expression.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_Traits.hpp"

#define UTOPIA_RECORDER_ENABLED
#ifdef UTOPIA_RECORDER_ENABLED
#define UTOPIA_RECORD_EXPR(macro_expr_, macro_tensor_) \
    { utopia::Recorder::instance().record_expr_and_value(macro_expr_, macro_tensor_); }
#define UTOPIA_RECORD_VALUE(macro_name_, macro_tensor_) \
    { utopia::Recorder::instance().record_name_and_value(macro_name_, macro_tensor_); }
#define UTOPIA_RECORD_SCOPE_BEGIN(macro_name_) \
    { utopia::Recorder::instance().scope_begin(macro_name_); }
#define UTOPIA_RECORD_SCOPE_END(macro_name_) \
    { utopia::Recorder::instance().scope_end(macro_name_); }
#else
#define UTOPIA_RECORD_EXPR(macro_expr_, macro_tensor_)
#define UTOPIA_RECORD_VALUE(macro_name_, macro_tensor_)
#define UTOPIA_RECORD_SCOPE_BEGIN(macro_name_)
#define UTOPIA_RECORD_SCOPE_END(macro_name_)
#endif  // UTOPIA_RECORDER_ENABLED

namespace utopia {
    class Recorder {
    public:
        template <class Derived, class T>
        void record_expr_and_value(const Expression<Derived> &expr, const Tensor<T, 1> &v) {
            record_name_and_value(expr.get_class(), v);
        }

        void scope_begin(const std::string &name) {
            ++n_nested_scopes_;

            if (mpi_world_rank() == 0) {
                os_ << std::setfill(' ');
                os_ << std::setw(n_nested_scopes_ * 10) << " "
                    << "%----------------------------------------------------------------\n";
                os_ << std::setw(n_nested_scopes_ * 10) << " "
                    << "%begin: " << name << std::endl;
            }
        }

        void scope_end(const std::string &name) {
            if (mpi_world_rank() == 0) {
                os_ << std::setfill(' ');
                os_ << std::setw(n_nested_scopes_ * 10) << " "
                    << "%end: " << name << std::endl;
                os_ << std::setw(n_nested_scopes_ * 10) << " "
                    << "%----------------------------------------------------------------\n";
            }
            --n_nested_scopes_;
        }

        template <class T>
        void record_name_and_value(const std::string &name, const Tensor<T, 1> &v) {
            typedef utopia::Tensor<T, 1> Vector;
            DEF_UTOPIA_SCALAR(Vector);

            mpi_world_barrier();

            if (mpi_world_rank() == 0) {
                os_ << std::setfill(' ');
                os_ << std::endl;
                os_ << std::setw(n_nested_scopes_ * 10) << " "
                    << "% " << name << "\n";
                os_ << std::setw(n_nested_scopes_ * 10) << " "
                    << "v_" << expr_num_++ << " = [";
            }

            mpi_world_barrier();

            for (int i = 0; i < mpi_world_size(); ++i) {
                mpi_world_barrier();

                if (i == mpi_world_rank()) {
                    os_ << std::flush;

                    each_read(v, [this](const SizeType /*i*/, const Scalar val) { os_ << val << " "; });

                    os_ << std::flush;
                }
            }

            mpi_world_barrier();

            os_ << std::flush;
            if (mpi_world_rank() == 0) {
                os_ << "]';";
                os_ << std::endl;
            }

            os_ << std::flush;

            mpi_world_barrier();
        }

        static Recorder &instance() {
            static Recorder instance_;
            return instance_;
        }

        Recorder(std::ostream &os = std::cout) : os_(os) {}

        std::ostream &os_;
        long expr_num_{0};
        long n_nested_scopes_{0};
    };
}  // namespace utopia

#endif  // UTOPIA_RECORDER_HPP

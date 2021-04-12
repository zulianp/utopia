#ifndef UTOPIA_STK_SPACE_IO_HPP
#define UTOPIA_STK_SPACE_IO_HPP

#include "utopia_Input.hpp"
#include "utopia_Path.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_stk_ForwardDeclarations.hpp"
#include "utopia_stk_FunctionSpace.hpp"

#include <string>

namespace utopia {
    namespace stk {

        class SpaceIO : public Configurable {
        public:
            using Scalar = Traits<FunctionSpace>::Scalar;
            using Vector = Traits<FunctionSpace>::Vector;

            void read(Input &in) override;
            bool read_with_state(Input &in, Field<FunctionSpace> &field);
            // bool load();

            bool write(const Vector &v);
            bool write(const Vector &v, const int step, const Scalar);
            // bool read(Vector &v, const int step = 1, const Scalar t = 0);

            void set_output_path(const Path &path);
            // void set_read_path(const Path &path);

            SpaceIO(FunctionSpace &space);
            ~SpaceIO();

        public:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace stk

    template <>
    class IO<utopia::stk::FunctionSpace> : public utopia::stk::SpaceIO {
    public:
        using utopia::stk::SpaceIO::SpaceIO;
    };

}  // namespace utopia

#endif  // UTOPIA_STK_SPACE_IO_HPP

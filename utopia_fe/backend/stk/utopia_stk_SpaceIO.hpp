#ifndef UTOPIA_STK_SPACE_IO_HPP
#define UTOPIA_STK_SPACE_IO_HPP

#include "utopia_Input.hpp"
#include "utopia_Path.hpp"

#include "utopia_Field.hpp"
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

            // One shot functions
            void read(Input &in) override;
            bool read_with_state(Input &in, Field<FunctionSpace> &field);
            bool read_nodal(Field<FunctionSpace> &field);
            void import_all_field_data(const bool value);

            // Statefull functions
            bool write(const Field<FunctionSpace> &field);
            bool write(const Vector &v);
            bool write(const Vector &v, const int step, const Scalar);
            bool write(const int step, const Scalar t);

            void set_output_path(const Path &path);
            bool open_output();

            void set_input_path(const Path &path);
            bool open_input();
            bool open_input(Input &in);

            SpaceIO(FunctionSpace &space);
            ~SpaceIO();

            void enable_interpolation_mode();

            bool load_time_step(const Scalar t);
            void register_output_field(const std::string &var_name);
            void register_output_field(const Field<FunctionSpace> &field);
            void update_output_field(const Field<FunctionSpace> &field);

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

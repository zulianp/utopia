#ifndef UTOPIA_PETSC_SPACE_IO_HPP
#define UTOPIA_PETSC_SPACE_IO_HPP

#include "utopia_Input.hpp"
#include "utopia_Path.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_FunctionSpace.hpp"

#include <string>

namespace utopia {
    namespace petsc {

        class SpaceIO : public Configurable {
        public:
            using Scalar = Traits<FunctionSpace>::Scalar;
            using Vector = Traits<FunctionSpace>::Vector;

            // One shot functions
            void read(Input &in) override;
            bool read_with_state(Input &in, Field<FunctionSpace> &field);
            void import_all_field_data(const bool value);

            // Statefull functions
            // bool load();

            bool write(const Vector &v);
            bool write(const Vector &v, const int step, const Scalar);
            bool write(const int step, const Scalar t);
            // bool read(Vector &v, const int step = 1, const Scalar t = 0);

            void set_output_path(const Path &path);
            bool open_output();

            void set_input_path(const Path &path);
            bool open_input();
            bool open_input(Input &in);
            // void set_read_path(const Path &path);

            SpaceIO(FunctionSpace &space);
            ~SpaceIO();

            void enable_interpolation_mode();

            bool load_time_step(const Scalar t);
            void register_output_field(const std::string &var_name);

        public:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace petsc

    template <>
    class IO<utopia::petsc::FunctionSpace> : public utopia::petsc::SpaceIO {
    public:
        using utopia::petsc::SpaceIO::SpaceIO;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_SPACE_IO_HPP

#ifndef UTOPIA_MARS_SPACE_IO_HPP
#define UTOPIA_MARS_SPACE_IO_HPP

#include "utopia_Input.hpp"
#include "utopia_Path.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_mars_ForwardDeclarations.hpp"
#include "utopia_mars_FunctionSpace.hpp"

#include <string>

namespace utopia {
    namespace mars {

        class SpaceIO : public Configurable {
        public:
            using Scalar = Traits<FunctionSpace>::Scalar;
            using Vector = Traits<FunctionSpace>::Vector;

            // One shot functions
            void read(Input &in) override { Utopia::Abort("IMPLEMENT ME!"); }
            bool read_with_state(Input &in, Field<FunctionSpace> &field) { Utopia::Abort("IMPLEMENT ME!"); }
            void import_all_field_data(const bool value) { Utopia::Abort("IMPLEMENT ME!"); }

            // Statefull functions
            // bool load() { Utopia::Abort("IMPLEMENT ME!");}

            bool write(const Vector &v) { Utopia::Abort("IMPLEMENT ME!"); }
            bool write(const Vector &v, const int step, const Scalar) { Utopia::Abort("IMPLEMENT ME!"); }
            bool write(const int step, const Scalar t) { Utopia::Abort("IMPLEMENT ME!"); }
            // bool read(Vector &v, const int step = 1, const Scalar t = 0) { Utopia::Abort("IMPLEMENT ME!");}

            void set_output_path(const Path &path) { Utopia::Abort("IMPLEMENT ME!"); }
            bool open_output() { Utopia::Abort("IMPLEMENT ME!"); }

            void set_input_path(const Path &path) { Utopia::Abort("IMPLEMENT ME!"); }
            bool open_input() { Utopia::Abort("IMPLEMENT ME!"); }
            bool open_input(Input &in) { Utopia::Abort("IMPLEMENT ME!"); }
            // void set_read_path(const Path &path) { Utopia::Abort("IMPLEMENT ME!");}

            SpaceIO(FunctionSpace &space) {}
            ~SpaceIO() {}

            void enable_interpolation_mode() { Utopia::Abort("IMPLEMENT ME!"); }

            bool load_time_step(const Scalar t) { Utopia::Abort("IMPLEMENT ME!"); }
            void register_output_field(const std::string &var_name) { Utopia::Abort("IMPLEMENT ME!"); }

        public:
            // class Impl;
            // std::unique_ptr<Impl> impl_;
        };

    }  // namespace mars

    template <>
    class IO<utopia::mars::FunctionSpace> : public utopia::mars::SpaceIO {
    public:
        using utopia::mars::SpaceIO::SpaceIO;
    };

}  // namespace utopia

#endif  // UTOPIA_MARS_SPACE_IO_HPP

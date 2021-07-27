#ifndef UTOPIA_LIBMESH_IO_HPP
#define UTOPIA_LIBMESH_IO_HPP

#include "utopia_Input.hpp"
#include "utopia_Path.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_libmesh_ForwardDeclarations.hpp"
#include "utopia_libmesh_FunctionSpace_new.hpp"

#include <string>

namespace utopia {
    namespace libmesh {

        class SpaceIO : public Configurable {
        public:
            using Scalar = Traits<FunctionSpace>::Scalar;
            using Vector = Traits<FunctionSpace>::Vector;

            void read(Input &in) override;
            bool read_with_state(Input &in, Field<FunctionSpace> &field);
            // bool load();

            bool write(const Vector &v);
            bool write(const Vector &v, const int step, const Scalar t);
            // bool read(Vector &v, const int step = 1, const Scalar t = 0);

            void set_output_path(const Path &path);
            // void set_read_path(const Path &path);

            SpaceIO(FunctionSpace &space);
            ~SpaceIO();

        public:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace libmesh

    template <>
    class IO<utopia::libmesh::FunctionSpace> : public utopia::libmesh::SpaceIO {
    public:
        using utopia::libmesh::SpaceIO::SpaceIO;
    };

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_IO_HPP

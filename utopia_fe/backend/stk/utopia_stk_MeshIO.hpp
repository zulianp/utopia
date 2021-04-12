#ifndef UTOPIA_STK_MESH_IO_HPP
#define UTOPIA_STK_MESH_IO_HPP

#include "utopia_Input.hpp"
#include "utopia_Path.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_stk_ForwardDeclarations.hpp"

#include <memory>
#include <string>

namespace utopia {
    namespace stk {

        class MeshIO : public Configurable {
        public:
            void read(Input &in) override;
            bool load();
            bool write(const Path &write_path);

            void set_read_path(const Path &path);
            void set_read_specification(const std::string &format);

            void import_all_field_data(const bool value);

            MeshIO(Mesh &mesh);
            ~MeshIO();

        public:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace stk

    template <>
    class IO<utopia::stk::Mesh> : public utopia::stk::MeshIO {
    public:
        using utopia::stk::MeshIO::MeshIO;
    };

}  // namespace utopia

#endif  // UTOPIA_STK_MESH_IO_HPP

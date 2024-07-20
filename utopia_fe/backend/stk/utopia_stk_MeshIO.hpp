#ifndef UTOPIA_STK_MESH_IO_HPP
#define UTOPIA_STK_MESH_IO_HPP

#include "utopia_Input.hpp"
#include "utopia_Path.hpp"
#include "utopia_Traits.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_stk_ForwardDeclarations.hpp"
#include "utopia_stk_Mesh.hpp"

#include <memory>
#include <string>

namespace utopia {

    namespace stk {

        class MeshIO : public Configurable {
        public:
            using Scalar = Traits<Mesh>::Scalar;

            void read(Input &in) override;
            bool load();
            bool load_time_step(const Scalar t);
            bool load_last_time_step();

            bool write(const Path &write_path);
            bool write(const int step, const Scalar t);

            void set_read_path(const Path &path);
            void set_read_specification(const std::string &format);

            void import_all_field_data(const bool value);

            ::stk::io::StkMeshIoBroker &raw_type();

            void set_output_path(const Path &path);
            void set_output_mode(enum OutputMode output_mode);
            void set_output_mode(const std::string &output_mode);

            void create_output_mesh();
            void register_output_field(const std::string &var_name);

            int num_time_steps() const;
            Scalar max_time() const;
            void set_import_all_data(const bool val);

            MeshIO(Mesh &mesh);
            ~MeshIO();

            const Path &output_path() const;
            bool ensure_output();
            bool ensure_input();

            void enable_interpolation_mode();

            std::vector<double> get_time_steps() const;

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

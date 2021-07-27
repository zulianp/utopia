
#ifndef UTOPIA_MARS_MESH_IO_HPP
#define UTOPIA_MARS_MESH_IO_HPP

// #include "utopia_Input.hpp"
// #include "utopia_Path.hpp"
// #include "utopia_Traits.hpp"

// #include "utopia_fe_Core.hpp"

// #include "utopia_mars_ForwardDeclarations.hpp"
// #include "utopia_mars_Mesh.hpp"

// #include <memory>
// #include <string>

#ifdef MARS_WITH_ADIOS2
#include "mars_adios2_IO.hpp"
#else
#ifdef MARS_WITH_VTK
#include "mars_vtk_IO.hpp"
#endif
#endif

namespace utopia {
    namespace mars {

#ifdef MARS_WITH_ADIOS2
#define MARS_WITH_WITH_IO
        template <class... Args>
        using MarsIOImpl = ::mars::adios2::IO<Args...>;
#else
#ifdef MARS_WITH_VTK
#define MARS_WITH_WITH_IO
        template <class... Args>
        using MarsIOImpl = ::mars::vtk::IO<Args...>;
#endif
#endif

        //         class MeshIO : public Configurable {
        //         public:
        //             using Scalar = Traits<Mesh>::Scalar;

        //             void read(Input &in) override;
        //             // bool load();
        //             // bool load_time_step(const Scalar t);
        //             bool write(const Path &write_path);
        //             // bool write(const int step, const Scalar t);

        //             // void set_read_path(const Path &path);
        //             // void set_read_specification(const std::string &format);

        //             // void import_all_field_data(const bool value);

        //             // void set_output_path(const Path &path);
        //             // void create_output_mesh();
        //             // void register_output_field(const std::string &var_name);

        //             MeshIO(Mesh &mesh);
        //             ~MeshIO();

        //             // const Path &output_path() const;
        //             // bool ensure_output();
        //             // bool ensure_input();

        //         public:
        //             class Impl;
        //             std::unique_ptr<Impl> impl_;
        //         };

        //     }  // namespace mars

        //     template <>
        //     class IO<utopia::mars::Mesh> : public utopia::mars::MeshIO {
        //     public:
        //         using utopia::mars::MeshIO::MeshIO;
        //     };
    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_MESH_IO_HPP

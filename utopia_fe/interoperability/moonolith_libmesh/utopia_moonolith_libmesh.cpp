#include "utopia_moonolith_libmesh.hpp"

// FIXME Legacy includes (they will have to be reorganized)
#include "utopia_libmesh_Mesh.hpp"
#include "utopia_moonolith_Mesh.hpp"

#include "utopia_libmesh_FunctionSpace_new.hpp"
#include "utopia_moonolith_FunctionSpace.hpp"

#include "utopia_moonolith_libmesh_Convert.hpp"

namespace utopia {

    void ConvertMesh<utopia::libmesh::Mesh, utopia::moonolith::Mesh>::apply(const utopia::libmesh::Mesh &in,
                                                                            utopia::moonolith::Mesh &out) {
        using Scalar = Traits<utopia::moonolith::Mesh>::Scalar;

        const int dim = in.spatial_dimension();

        switch (dim) {
            case 1: {
                auto m_mesh = std::make_shared<::moonolith::Mesh<Scalar, 1>>(in.comm().raw_comm());
                convert(in.raw_type(), *m_mesh);
                out.wrap(m_mesh);
                break;
            }

            case 2: {
                auto m_mesh = std::make_shared<::moonolith::Mesh<Scalar, 2>>(in.comm().raw_comm());
                convert(in.raw_type(), *m_mesh);
                out.wrap(m_mesh);
                break;
            }

            case 3: {
                auto m_mesh = std::make_shared<::moonolith::Mesh<Scalar, 3>>(in.comm().raw_comm());
                convert(in.raw_type(), *m_mesh);
                out.wrap(m_mesh);
                break;
            }
            default: {
                Utopia::Abort();
            }
        }
    }

    void ConvertFunctionSpace<utopia::libmesh::FunctionSpace, utopia::moonolith::FunctionSpace>::apply(
        const utopia::libmesh::FunctionSpace &in,
        utopia::moonolith::FunctionSpace &out) {
        using Scalar = Traits<utopia::moonolith::Mesh>::Scalar;

        const int dim = in.mesh().spatial_dimension();

        switch (dim) {
            case 1: {
                auto m_mesh = std::make_shared<::moonolith::Mesh<Scalar, 1>>(in.comm().raw_comm());
                auto m_space = std::make_shared<::moonolith::FunctionSpace<::moonolith::Mesh<Scalar, 1>>>(m_mesh);
                convert_libmesh_to_moonolith(in.raw_type().get_mesh(), in.raw_type_dof_map(), 0, *m_space, 0);
                out.wrap(m_space);
                break;
            }

            case 2: {
                auto m_mesh = std::make_shared<::moonolith::Mesh<Scalar, 2>>(in.comm().raw_comm());
                auto m_space = std::make_shared<::moonolith::FunctionSpace<::moonolith::Mesh<Scalar, 2>>>(m_mesh);
                convert_libmesh_to_moonolith(in.raw_type().get_mesh(), in.raw_type_dof_map(), 0, *m_space, 0);
                out.wrap(m_space);
                break;
            }

            case 3: {
                auto m_mesh = std::make_shared<::moonolith::Mesh<Scalar, 3>>(in.comm().raw_comm());
                auto m_space = std::make_shared<::moonolith::FunctionSpace<::moonolith::Mesh<Scalar, 3>>>(m_mesh);
                convert_libmesh_to_moonolith(in.raw_type().get_mesh(), in.raw_type_dof_map(), 0, *m_space, 0);
                out.wrap(m_space);
                break;
            }

            default: {
                Utopia::Abort();
            }
        }
    }

    void ExtractTraceSpace<utopia::libmesh::FunctionSpace, utopia::moonolith::FunctionSpace>::apply(
        const utopia::libmesh::FunctionSpace &in,
        utopia::moonolith::FunctionSpace &out) {
        using Scalar = Traits<utopia::libmesh::FunctionSpace>::Scalar;

        const int dim = in.mesh().spatial_dimension();

        switch (in.mesh().spatial_dimension()) {
            case 2: {
                using Mesh2_t = ::moonolith::Mesh<Scalar, 2>;
                using Space2_t = ::moonolith::FunctionSpace<Mesh2_t>;

                auto mesh = std::make_shared<Mesh2_t>(in.comm().raw_comm());
                auto temp_space = std::make_shared<Space2_t>(mesh);

                extract_trace_space(in.mesh().raw_type(), in.raw_type_dof_map(), 0, *temp_space, {});

                out.wrap(temp_space);
                break;
            }

            case 3: {
                using Mesh3_t = ::moonolith::Mesh<Scalar, 3>;
                using Space3_t = ::moonolith::FunctionSpace<Mesh3_t>;

                auto mesh = std::make_shared<Mesh3_t>(in.comm().raw_comm());
                auto temp_space = std::make_shared<Space3_t>(mesh);

                extract_trace_space(in.mesh().raw_type(), in.raw_type_dof_map(), 0, *temp_space, {});

                out.wrap(temp_space);
                break;
            }

            default: {
                assert(false);
                Utopia::Abort("ExtractTraceSpace<FunctionSpace> spatial_dimension not supported");
            }
        }
    }

}  // namespace utopia

#ifndef UTOPIA_LIBMESH_DOF_MAP_ADAPTER_HPP
#define UTOPIA_LIBMESH_DOF_MAP_ADAPTER_HPP

#include "utopia_ElementDofMap.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "moonolith_communicator.hpp"

#include <vector>
#include <memory>

#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/boundary_mesh.h"
#include "libmesh/dof_map.h"

namespace utopia {

    class LibMeshDofMapAdapter {
    public:

        inline void describe(std::ostream &os = std::cout) const
        {
            os << "n_ne_dofs: " << mapping_.from_extent() << "\n";
            os << "element_dof_map: [" << mapping_.from_range_begin << "-" << mapping_.from_range_end << ")\n";
            os << "dof_map:         [" << mapping_.to_range_begin   << "-" << mapping_.to_range_end   << ")\n";

            if(permutation_) {
                auto perm_size = size(*permutation_);
                double entries = sum(*permutation_);

                std::cout << perm_size.get(0) << " x "  << perm_size.get(1) << " nnz: " << entries << std::endl;
            }

            if(vector_permutation_) {
                auto perm_size = size(*vector_permutation_);
                double entries = sum(*vector_permutation_);

                std::cout << perm_size.get(0) << " x "  << perm_size.get(1) << " nnz: " << entries << std::endl;
            }

        }

        class Mapping {
        public:
            SizeType from_range_begin, from_range_end;
            SizeType to_range_begin, to_range_end;
            std::vector<SizeType> idx;

            inline std::size_t size() const { return idx.size(); }
            inline SizeType from_extent() const { return from_range_end - from_range_begin; }
            inline SizeType to_extent()   const { return to_range_end   - to_range_begin; }
        };
        
        Mapping & mapping()
        {
            return mapping_;
        }
        
        std::vector<ElementDofMap> & element_dof_map()
        {
            return dofs_;
        }

        const std::vector<ElementDofMap> & element_dof_map() const
        {
            return dofs_;
        }
        
        std::shared_ptr<USparseMatrix> & permutation()
        {
            return permutation_;
        }

        const Mapping & mapping() const
        {
            return mapping_;
        }
        
        const std::vector<ElementDofMap> & dofs() const
        {
            return dofs_;
        }
        
        std::shared_ptr<const USparseMatrix> permutation() const
        {
            assert(permutation_);
            return permutation_;
        }

        std::shared_ptr<const USparseMatrix> vector_permutation() const
        {
            assert(vector_permutation_);
            return vector_permutation_;
        }

        inline SizeType n_local_dofs() const
        {
            return mapping_.from_extent();
        }

        //surf_mesh is extracted from vol_mesh
        void init_for_contact(
            libMesh::MeshBase &vol_mesh,
            const libMesh::DofMap &vol_dof_map,
            libMesh::MeshBase &surf_mesh,
            const libMesh::DofMap &surf_dof_map,
            const int var_num,
            const int surf_var_num);

        void init(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const libMesh::DofMap &dof_map,
            const int var_num,
            const bool build_element_dof_map);

        //let us assume that the displacement degees of freedom are consecutive 
        //in the volume dofmap
        static void vector_permuation_map_from_map(
            const std::size_t spatial_dim,
            const Mapping &mapping,
            Mapping &vector_mapping);

        static void permutation_matrix_from_map(
            const Mapping &map,
            USparseMatrix &mat);

        static void boundary_permutation_map(
            const libMesh::MeshBase &surf_mesh,
            const libMesh::DofMap   &volume_dof_map,
            const libMesh::DofMap   &surf_dof_map,
            unsigned int var_num,
            unsigned int surf_var_num,
            Mapping &map);

        static void boundary_element_node_permutation_map(
            const libMesh::MeshBase &surf_mesh,
            const libMesh::DofMap   &volume_dof_map,
            const libMesh::DofMap   &surf_dof_map,
            unsigned int var_num,
            unsigned int surf_var_num,
            Mapping &map);

        static SizeType count_dof_x_elem(
            const libMesh::MeshBase &mesh,
            const libMesh::DofMap &dof_map
            );

        static Range element_node_range(
            const libMesh::MeshBase &mesh,
            const libMesh::DofMap &dof_map);

        static void element_node_dof_map_and_permutation(
            const libMesh::MeshBase &mesh,
            const libMesh::DofMap &dof_map,
            const Mapping &map,
            std::vector<ElementDofMap> &elem_dof_map,
            USparseMatrix &mat);

    private:
        Mapping mapping_;
        std::vector<ElementDofMap> dofs_;
        std::shared_ptr<USparseMatrix> permutation_;

        //only for contact
        Mapping vector_mapping_;
        std::shared_ptr<USparseMatrix> vector_permutation_;
    };
}

#endif //UTOPIA_LIBMESH_DOF_MAP_ADAPTER_HPP

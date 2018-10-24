#ifndef UTOPIA_TRANSFER_LOCAL_2_GLOBAL_ASSEMBLER_HPP
#define UTOPIA_TRANSFER_LOCAL_2_GLOBAL_ASSEMBLER_HPP

#include "utopia.hpp"
#include "utopia_Grid.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_QMortarBuilder.hpp"
#include "MortarAssemble.hpp"

#include "utopia_intersector.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_Socket.hpp"

#include "libmesh/serial_mesh.h"
#include <functional>
#include <array>
#include <vector>

namespace utopia {
    namespace private_ {
        
        class PetrovGalerkinAssembler {
        public:
            using Elem   = libMesh::Elem;
            using FEType = libMesh::FEType;
            using Matrix = libMesh::DenseMatrix<libMesh::Real>;

            using FunctionSpace = utopia::LibMeshFunctionSpace;
            using SparseMatrix  = utopia::USparseMatrix;
            using MeshBase      = libMesh::MeshBase;
            using DofMap        = libMesh::DofMap;

            using ElementMatrix = LocalAssembler::Matrix;
            
            PetrovGalerkinAssembler();
            
            //dof_fun(master_dofs, slave_dofs)
            bool assemble(const Elem &master,
                          FEType master_type,
                          const Elem &slave,
                          FEType slave_type,
                          std::function<void(std::vector<long> &, std::vector<long> &)> dof_fun
            );
            
            void initialize(const moonolith::Communicator &comm,
                            const std::shared_ptr<LocalAssembler> &assembler,
                            const std::shared_ptr<Local2Global> &local2global,
                            const TransferOptions &opts,
                            const SizeType from_n_dofs,
                            const SizeType from_n_local_dofs,
                            const SizeType to_n_dofs,
                            const SizeType to_n_local_dofs);
            
            
            void finalize(std::vector<std::shared_ptr<SparseMatrix>> &mats);
            void print_stats();
            
        private:
            
            moonolith::Communicator comm_;
            std::shared_ptr<LocalAssembler> assembler_;
            std::shared_ptr<Local2Global> local2global_;

            TransferOptions opts_;

            SizeType from_n_dofs_;
            SizeType from_n_local_dofs_;
            SizeType to_n_dofs_;
            SizeType to_n_local_dofs_;
            long n_intersections_;
            long n_false_positives_;
            
            std::vector< libMesh::DenseMatrix<libMesh::Real> > elemmat_;
            std::vector< libMesh::Real > local_element_matrices_sum_;
            std::vector< std::shared_ptr< moonolith::SparseMatrix<double> > > mat_buffer_;
            
            std::vector<long> master_dofs_, slave_dofs_;
            
            void init_buffers();
            void finalize_form(std::size_t buffer_num, SparseMatrix &mat);
        };

    }
}

#endif //UTOPIA_TRANSFER_LOCAL_2_GLOBAL_ASSEMBLER_HPP

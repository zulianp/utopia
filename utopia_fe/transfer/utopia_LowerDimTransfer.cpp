#include "utopia_LowerDimTransfer.hpp"


#include "utopia_MatrixInserter.hpp"
#include "utopia_ZeroRowsToIdentity.hpp"
#include "utopia_NormalizeRows.hpp"
#include "utopia_LibMeshShape.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"


#include "moonolith_affine_transform.hpp"
#include "moonolith_contact.hpp"
#include "moonolith_sparse_matrix.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_assign_functions.hpp"
#include "moonolith_elem_shape.hpp"
#include "moonolith_elem_triangle.hpp"
#include "moonolith_elem_quad.hpp"
#include "moonolith_elem_segment.hpp"
#include "moonolith_matlab_scripter.hpp"

#include "par_moonolith_config.hpp"
#include "moonolith_sparse_matrix.hpp"
#include "moonolith_duplicate_intersection_avoidance.hpp"
#include "moonolith_one_master_one_slave_algorithm.hpp"

#include "moonolith_l2_assembler.hpp"
#include "moonolith_keast_quadrature_rule.hpp"
#include "moonolith_gauss_quadrature_rule.hpp"
#include "moonolith_par_l2_transfer.hpp"

#include "moonolith_par_lower_dim_l2_transfer.hpp"

#include "utopia_LibMeshFunctionSpaceAdapter.hpp"
#include "utopia_LibMeshToMoonolithConvertions.hpp"
#include "utopia_moonolith_permutations.hpp"
#include "utopia_NormalizeRows.hpp"
#include "utopia_TransferUtils.hpp"

namespace utopia {

    void LowerDimTransfer::TransferData::permute(const USparseMatrix &P, TransferData &out)
    {
        if(!empty(*B)) *out.B = P * *B;
        if(!empty(*D)) *out.D = P * *D * transpose(P);
        
        if(!empty(*T)) *out.T = P * *T;

        // if(!empty(*Q)) {
        //     *out.Q = P * *Q * transpose(P);
        //     normalize_rows(*out.Q);
        // }
    }

    using TransferDataT = utopia::LowerDimTransfer::TransferData;

    template<int Dim>
    class LowerDimTransferAlgorithm {
    public:
        using MeshT = moonolith::Mesh<double, Dim>;
        using FunctionSpaceT = moonolith::FunctionSpace<MeshT>;

        template<class Transfer>
        static void prepare_data(
            const TransferOptions &opts,
            Transfer &t,
            TransferDataT &data)
        {
            auto &B = *data.B;
            auto &D = *data.D;
            // auto &Q = *data.Q;
            auto &T = *data.T;

            convert_matrix(t.buffers.B.get(), B);
            convert_matrix(t.buffers.D.get(), D);
            // convert_matrix(t.buffers.Q.get(), Q);

            // if(!empty(Q)) {
                // m_utopia_warning_once("using sum(D, 1) instead of diag(D)");

                //pseudo-stuff
                UVector d_inv = sum(D, 1);

                e_pseudo_inv(d_inv, d_inv, 1e-12);

                USparseMatrix D_tilde_inv = diag(d_inv);
                USparseMatrix T_temp = D_tilde_inv * B;

                if(opts.n_var == 1) {
                    // T = Q * T_temp;
                    T =  T_temp;
                } else {
                    // USparseMatrix T_x = Q * T_temp;
                    // tensorize(T_x, opts.n_var, T);

                    tensorize(T_temp, opts.n_var, T);
                }
            // }
        }

        static bool apply(
            const libMesh::MeshBase &lm_mesh,
            const libMesh::DofMap   &lm_dofs,
            const TransferOptions &opts,
            TransferDataT &data)
        {
            moonolith::Communicator comm = lm_mesh.comm().get();
            auto mesh = std::make_shared<MeshT>(comm);

            FunctionSpaceT space(mesh);
            convert(lm_mesh, lm_dofs, opts.from_var_num, space);

            comm.barrier();

            if(comm.is_root()) {
                moonolith::logger() << "LowerDimTransferAlgorithm:surface_apply(...) begin" << std::endl;
            }
            
            static const int ManifoldDim = moonolith::StaticMax<Dim-1, 1>::value;

            moonolith::ParLowerDimL2Transfer<
                double,
                Dim,
                ManifoldDim,
                ManifoldDim
            > assembler(comm);

            if(!assembler.assemble(space, opts.tags, 1e-8)) {
                return false;
            }

            prepare_data(opts, assembler, data);

            comm.barrier();
            
            if(comm.is_root()) {
                moonolith::logger() << "LowerDimTransferAlgorithm:surface_apply(...) end" << std::endl;
            }

            return true;
        }
    };

    bool LowerDimTransfer::assemble(
        const MeshBase &mesh,
        const DofMap   &dofs,
        const TransferOptions &opts
    )
    {
        auto spatial_dim = mesh.spatial_dimension();
        bool has_intersection = false;

        if(spatial_dim == 2) {
            has_intersection = LowerDimTransferAlgorithm<2>::apply(mesh, dofs, opts, data);
        } else if(spatial_dim == 3) {
            has_intersection = LowerDimTransferAlgorithm<3>::apply(mesh, dofs, opts, data);
        } else {
            assert(false && "cannot be used with this dimension");
        }

        return has_intersection;
    }
}


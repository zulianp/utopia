#include "utopia_NewTransferAssembler.hpp"

#include "utopia_DualBasis.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_LibMeshBackend.hpp"
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


#include "utopia_LibMeshFunctionSpaceAdapter.hpp"


namespace utopia {

    template<int Dim>
    class TransferAlgorithm {
    public:
        using MasterAndSlaveAlgorithmT = moonolith::OneMasterOneSlaveAlgorithm<Dim, LibMeshFunctionSpaceAdapter>;
        using Adapter                  = typename MasterAndSlaveAlgorithmT::Adapter;
        // using SingleSpaceAlgorithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, SpaceT>;

        using Elem = moonolith::Elem<double, Dim, Dim>;
        using QuadratureT = moonolith::Quadrature<double, Dim>;
        using L2TransferT = moonolith::L2Transfer<Elem, Elem>;

        class LocalAssembler {
        public:

            inline bool operator()(const Adapter &master, const Adapter &slave)
            {
                auto &e_m = master.elem();
                auto &e_s = slave.elem();

                auto &m_m = master.collection();
                auto &m_s = slave.collection();

                auto &dofs_m = master.dofs().dofs;
                auto &dofs_s = slave.dofs().dofs;

                make(e_m, m_m.libmesh_fe_type(0), master_elem);
                make(e_s, m_s.libmesh_fe_type(0), slave_elem);

                init_quadrature(
                    order_for_l2_integral(
                        Dim,
                        e_m,
                        m_m.fe_type(0).order,
                        e_s,
                        m_s.fe_type(0).order
                    )
                );

                if(algo.assemble(*master_elem, *slave_elem)) {
                    uint n_nodes_master = dofs_m.size();
                    uint n_nodes_slave  = dofs_s.size();
                    
                    const auto &B_e = algo.coupling_matrix();
                    const auto &D_e = algo.mass_matrix();
                    const auto &Q_e = algo.transformation();

                    // std::cout << "---------------------\n";
                    // D_e.describe(std::cout);
                    // std::cout << "---------------------\n";

                    B.insert(dofs_s, dofs_m, B_e);
                    D.insert(dofs_s, dofs_s, D_e);
                    Q.insert(dofs_s, dofs_s, Q_e);

                    return true;
                } else {
                    return false;
                }
            }

            void init_quadrature(const int order)
            {
                moonolith::Gauss::get(order, q_rule);
                algo.set_quadrature(q_rule);
            }

            LocalAssembler(moonolith::Communicator &comm)
            : B(comm.get()), D(comm.get()), Q(comm.get(), false)
            {}

            QuadratureT q_rule;
            L2TransferT algo;
            std::shared_ptr<Elem> master_elem, slave_elem;

            MatrixInserter B, D, Q;
        };

        static bool apply(
            LibMeshFunctionSpaceAdapter &master,
            LibMeshFunctionSpaceAdapter &slave,
            const TransferOptions &opts,
            TransferData &data)
        {
            Chrono c;
            c.start();

            moonolith::Communicator comm(master.comm().get());

            MasterAndSlaveAlgorithmT algo(comm,
                moonolith::make_unique<LibMeshCollectionManagerT>(master.comm(), nullptr, true)
            );
            
            algo.init_simple(
                master,
                slave,
                0.0
            );

            c.stop();
            logger() << "init: " << c << std::endl;

            ////////////////////////////////////////////////////
            /////////////////// pair-wise method ///////////////

            c.start();
            
            LocalAssembler assembler(comm);
            algo.compute([&](const Adapter &master, const Adapter &slave) -> bool {
                if(assembler(master, slave)) {
                    return true;
                }

                return false;
            });

            double vol = assembler.algo.intersection_measure();

            comm.all_reduce(&vol, 1, moonolith::MPISum());

            double sum_mat_B = assembler.B.m_matrix.sum();
            double sum_mat_D = assembler.D.m_matrix.sum();

            auto n_master_dofs = master.n_local_dofs();
            auto n_slave_dofs  = slave.n_local_dofs();

            assembler.B.finalize(n_slave_dofs, n_master_dofs);
            assembler.D.finalize(n_slave_dofs, n_slave_dofs);
            assembler.Q.finalize(n_slave_dofs, n_slave_dofs);

            auto &B = *data.B;
            auto &D = *data.D;
            auto &Q = *data.Q;
            auto &T = *data.T;

            assembler.B.fill(B);
            assembler.D.fill(D);
            assembler.Q.fill(Q);

            if(!empty(Q)) {
                UVector d_inv = diag(D);
                e_pseudo_inv(d_inv, d_inv, 1e-15);

                USparseMatrix D_tilde_inv = diag(d_inv);
                USparseMatrix D_inv = Q * D_tilde_inv;
                // T = D_inv * B;
                T = D_tilde_inv * B;

                B.implementation().set_name("b");
                D.implementation().set_name("d");
                Q.implementation().set_name("q");
                T.implementation().set_name("t");

                write("B.m", B);
                write("D.m", D);
                write("Q.m", Q);
                write("T.m", T);

                // normalize_rows(T);
            }

            c.stop();
            logger() << "time MarsMeshTransfer::assemble: " << c  << std::endl;
            logger() << "vol: " << vol << " sum(B): " << sum_mat_B << " sum(D): " << sum_mat_D << std::endl;
            return vol > 0.0;
        }
    };


    bool NewTransferAssembler::assemble(
        const std::shared_ptr<MeshBase> &from_mesh,
        const std::shared_ptr<DofMap>   &from_dofs,
        const std::shared_ptr<MeshBase> &to_mesh,
        const std::shared_ptr<DofMap>   &to_dofs,
        const TransferOptions &opts
    )
    {
        LibMeshFunctionSpaceAdapter from_adapter, to_adapter;
        from_adapter.init(from_mesh, *from_dofs, opts.from_var_num);
        to_adapter.init(to_mesh, *to_dofs, opts.to_var_num);

        auto spatial_dim = to_mesh->spatial_dimension();

        bool has_intersection = false;

        if(spatial_dim == 1) {
            has_intersection = TransferAlgorithm<1>::apply(from_adapter, to_adapter, opts, data);
        } else if(spatial_dim == 2) {
            has_intersection = TransferAlgorithm<2>::apply(from_adapter, to_adapter, opts, data);
        } else if(spatial_dim == 3) {
            has_intersection = TransferAlgorithm<3>::apply(from_adapter, to_adapter, opts, data);
        }

        return has_intersection;
    }
}


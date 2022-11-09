#include "utopia_NewTransferAssembler.hpp"

#include "utopia_DualBasis.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_LibMeshShape.hpp"
#include "utopia_MatrixInserter.hpp"
#include "utopia_NormalizeRows.hpp"
#include "utopia_ZeroRowsToIdentity.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"

#include "moonolith_affine_transform.hpp"
#include "moonolith_assign_functions.hpp"
#include "moonolith_contact.hpp"
#include "moonolith_elem_quad.hpp"
#include "moonolith_elem_segment.hpp"
#include "moonolith_elem_shape.hpp"
#include "moonolith_elem_triangle.hpp"
#include "moonolith_matlab_scripter.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_sparse_matrix.hpp"

#include "moonolith_duplicate_intersection_avoidance.hpp"
#include "moonolith_one_master_one_slave_algorithm.hpp"
#include "moonolith_sparse_matrix.hpp"
#include "par_moonolith_config.hpp"

#include "moonolith_gauss_quadrature_rule.hpp"
#include "moonolith_keast_quadrature_rule.hpp"
#include "moonolith_l2_assembler.hpp"
#include "moonolith_par_l2_assembler_no_covering.hpp"
#include "moonolith_par_l2_transfer.hpp"
#include "moonolith_par_volume_surface_l2_transfer.hpp"

#include "utopia_Adaptivity.hpp"
#include "utopia_LibMeshFunctionSpaceAdapter.hpp"
#include "utopia_NormalizeRows.hpp"
#include "utopia_TransferUtils.hpp"
#include "utopia_moonolith_libmesh_Convert.hpp"
#include "utopia_moonolith_permutations.hpp"

namespace utopia {

    void NewTransferAssembler::TransferData::permute(const USparseMatrix &P, TransferData &out) {
        if (!empty(*B)) *out.B = P * *B;
        if (!empty(*D)) *out.D = P * *D * transpose(P);

        if (!empty(*T)) *out.T = P * *T;

        if (!empty(*Q)) {
            *out.Q = P * *Q * transpose(P);
            normalize_rows(*out.Q);
        }
    }

    using TransferDataT = utopia::NewTransferAssembler::TransferData;

    template <int Dim>
    class TransferAlgorithm {
    public:
        using MasterAndSlaveAlgorithmT = moonolith::OneMasterOneSlaveAlgorithm<Dim, LibMeshFunctionSpaceAdapter>;
        using Adapter = typename MasterAndSlaveAlgorithmT::Adapter;
        // using SingleSpaceAlgorithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, SpaceT>;

        using Elem = moonolith::Elem<double, Dim, Dim>;
        using QuadratureT = moonolith::Quadrature<double, Dim>;
        using L2TransferT = moonolith::L2Transfer<Elem, Elem>;

        class LocalAssembler {
        public:
            inline bool operator()(const Adapter &master, const Adapter &slave) {
                auto &e_m = master.elem();
                auto &e_s = slave.elem();

                auto &m_m = master.collection();
                auto &m_s = slave.collection();

                auto &dofs_m = master.dofs().dofs;
                auto &dofs_s = slave.dofs().dofs;

                make(e_m, m_m.libmesh_fe_type(0), master_elem);
                make(e_s, m_s.libmesh_fe_type(0), slave_elem);

                init_quadrature(order_for_l2_integral(Dim, e_m, m_m.fe_type(0).order, e_s, m_s.fe_type(0).order));

                if (algo.assemble(*master_elem, *slave_elem)) {
                    uint n_nodes_master = dofs_m.size();
                    uint n_nodes_slave = dofs_s.size();

                    const auto &B_e = algo.coupling_matrix();
                    const auto &D_e = algo.mass_matrix();
                    const auto &Q_e = algo.transformation();

                    B.insert(dofs_s, dofs_m, B_e);
                    D.insert(dofs_s, dofs_s, D_e);
                    Q.insert(dofs_s, dofs_s, Q_e);

                    return true;
                } else {
                    return false;
                }
            }

            void init_quadrature(const int order) {
                moonolith::Gauss::get(order, q_rule);
                algo.set_quadrature(q_rule);
            }

            LocalAssembler(moonolith::Communicator &comm) : B(comm.get()), D(comm.get()), Q(comm.get(), false) {}

            QuadratureT q_rule;
            L2TransferT algo;
            std::shared_ptr<Elem> master_elem, slave_elem;

            MatrixInserter B, D, Q;
        };

        static bool apply(LibMeshFunctionSpaceAdapter &master,
                          LibMeshFunctionSpaceAdapter &slave,
                          const TransferOptions &opts,
                          TransferDataT &data) {
            Chrono c;
            c.start();

            moonolith::Communicator comm(master.comm().get());
            moonolith::SearchSettings settings;

            if (Utopia::instance().verbose()) {
                // moonolith::root_describe("---------------------------------------\n"
                //     "begin: search_and_compute ",
                //     comm, std::cout);

                // settings.verbosity_level = 2;
                // settings.disable_redistribution = true;
            }

            double vol = 0.0;
            LocalAssembler assembler(comm);

            {
                MasterAndSlaveAlgorithmT algo(
                    comm, moonolith::make_unique<LibMeshCollectionManagerT>(master.comm(), nullptr, true), settings);

                if (opts.tags.empty()) {
                    algo.init_simple(master, slave, 0.0);
                } else {
                    algo.init(master, slave, opts.tags, 0.0);
                }

                c.stop();
                logger() << "init: " << c << std::endl;

                ////////////////////////////////////////////////////
                /////////////////// pair-wise method ///////////////

                c.start();

                algo.compute([&](const Adapter &master, const Adapter &slave) -> bool {
                    if (assembler(master, slave)) {
                        return true;
                    }

                    return false;
                });

                vol = assembler.algo.intersection_measure();
            }

            comm.all_reduce(&vol, 1, moonolith::MPISum());

            double sum_mat_B = assembler.B.m_matrix.sum();
            double sum_mat_D = assembler.D.m_matrix.sum();

            auto n_master_dofs = master.n_local_dofs();
            auto n_slave_dofs = slave.n_local_dofs();

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

            if (!empty(Q)) {
                UVector d_inv = sum(D, 1);

                // write("d_inv.m", d_inv);

                e_pseudo_inv(d_inv, d_inv, 1e-12);

                USparseMatrix D_tilde_inv = diag(d_inv);
                USparseMatrix T_temp = D_tilde_inv * B;

                if (opts.n_var == 1) {
                    T = Q * T_temp;
                } else {
                    USparseMatrix T_x = Q * T_temp;
                    tensorize(T_x, opts.n_var, T);
                }

                // B.implementation().set_name("b");
                // D.implementation().set_name("d");
                // Q.implementation().set_name("q");
                // T.implementation().set_name("t");

                // write("B.m", B);
                // write("D.m", D);
                // write("Q.m", Q);
                // write("T.m", T);

                // normalize_rows(T);
            }

            c.stop();
            logger() << "time MarsMeshTransfer::assemble: " << c << std::endl;
            logger() << "vol: " << vol << " sum(B): " << sum_mat_B << " sum(D): " << sum_mat_D << std::endl;
            return vol > 0.0;
        }
    };

    static void handle_constraints_pre_process(TransferDataT &data) {
        assert(data.constraint_matrix_from);
        assert(data.constraint_matrix_to);
        assert(data.D);
        assert(data.B);

        auto &cm_from = *data.constraint_matrix_from;
        auto &cm_to = *data.constraint_matrix_to;
        auto &D = *data.D;
        auto &B = *data.B;

        // rename("b", B);
        // write("B.m", B);

        // rename("d", D);
        // write("D.m", D);

        if (!empty(cm_to)) {
            // rename("cm_to", cm_to);
            // write("CM_to.m", cm_to);
            D = transpose(cm_to) * D * (cm_to);

            if (!empty(cm_from)) {
                // rename("cm_from", cm_from);
                // write("CM_from.m", cm_from);
                B = transpose(cm_to) * B * (cm_from);
            } else {
                B = transpose(cm_to) * B;
            }

        } else {
            if (!empty(cm_from)) {
                // rename("b", B);
                // write("B.m", B);

                // rename("cm", cm_from);
                // write("CM.m", cm_from);

                B = B * (cm_from);
            }
        }
    }

    static void handle_constraints_post_process(TransferDataT &data, USparseMatrix &temp_T) {
        auto &post_constraint_matrix_to = *data.post_constraint_matrix_to;

        if (empty(post_constraint_matrix_to)) return;

        temp_T += post_constraint_matrix_to * temp_T;

        // rename("tc", temp_T);
        // write("TC.m", temp_T);
        //
        //        rename("hn", post_constraint_matrix_to);
        //        write("HN.m", post_constraint_matrix_to);
    }

    template <class Transfer>
    static void prepare_data(const TransferOptions &opts, Transfer &t, TransferDataT &data) {
        auto &B = *data.B;
        auto &D = *data.D;
        auto &Q = *data.Q;
        auto &T = *data.T;

        convert_matrix(t.buffers.B.get(), B);
        convert_matrix(t.buffers.D.get(), D);
        convert_matrix(t.buffers.Q.get(), Q);

        if (!empty(Q)) {
            m_utopia_warning_once("using sum(D, 1) instead of diag(D)");

            handle_constraints_pre_process(data);

            UVector d_inv = sum(D, 1);

            e_pseudo_inv(d_inv, d_inv, 1e-12);

            USparseMatrix D_tilde_inv = diag(d_inv);
            USparseMatrix T_x = Q * D_tilde_inv * B;

            handle_constraints_post_process(data, T_x);

            if (opts.n_var == 1) {
                T = T_x;
            } else {
                tensorize(T_x, opts.n_var, T);
            }

            // rename("t", T);
            // write("T.m", T);
        }
    }

    template <int Dim, int DimMaster, int DimSlave>
    class ApplyAux final {
    public:
        using MeshT = moonolith::Mesh<double, Dim>;
        using FunctionSpaceT = moonolith::FunctionSpace<MeshT>;

        static bool apply(const FunctionSpaceT &master,
                          const FunctionSpaceT &slave,
                          const TransferOptions &opts,
                          TransferDataT &data) {
            moonolith::Communicator comm = master.mesh().comm();
            comm.barrier();

            if (comm.is_root()) {
                moonolith::logger() << "ConvertTransferAlgorithm:apply(...) begin" << std::endl;
            }

            if (opts.has_covering) {
                moonolith::ParL2Transfer<double,
                                         Dim,
                                         moonolith::StaticMax<DimMaster, 1>::value,
                                         moonolith::StaticMax<DimSlave, 1>::value>
                    assembler(comm);

                if (opts.tags.empty()) {
                    if (!assembler.assemble(master, slave)) {
                        return false;
                    }
                } else {
                    if (!assembler.assemble(master, slave, opts.tags)) {
                        return false;
                    }
                }

                prepare_data(opts, assembler, data);

            } else {
                moonolith::ParL2TransferNoCovering<double,
                                                   Dim,
                                                   moonolith::StaticMax<DimMaster, 1>::value,
                                                   moonolith::StaticMax<DimSlave, 1>::value>
                    assembler(comm);

                // if (opts.tags.empty()) {
                if (!assembler.assemble(master, slave)) {
                    return false;
                }
                // } else {
                //     if (!assembler.assemble(master, slave, opts.tags)) {
                //         return false;
                //     }
                // }

                prepare_data(opts, assembler, data);
            }

            comm.barrier();

            if (comm.is_root()) {
                moonolith::logger() << "ConvertTransferAlgorithm:apply(...) end" << std::endl;
            }

            return true;
        }
    };

    template <int Dim>
    class ApplyAuxLowerDim {
    public:
        using MeshT = moonolith::Mesh<double, Dim>;
        using FunctionSpaceT = moonolith::FunctionSpace<MeshT>;

        static bool apply(const FunctionSpaceT &master,
                          const FunctionSpaceT &slave,
                          const TransferOptions &opts,
                          TransferDataT &data) {
            moonolith::Communicator comm = master.mesh().comm();
            comm.barrier();

            if (comm.is_root()) {
                moonolith::logger() << "ConvertTransferAlgorithm::apply(...)/Vol2Surf begin" << std::endl;
            }

            moonolith::ParVolumeSurfaceL2Transfer<double, Dim> assembler(comm);

            if (opts.tags.empty()) {
                if (!assembler.assemble(master, slave)) {
                    return false;
                }
            } else {
                // if(!assembler.assemble(master, slave, opts.tags)) {
                //     return false;
                // }
                std::cerr << "[Error] not implemented" << std::endl;
                assert(false && "implement me!!!");
            }

            prepare_data(opts, assembler, data);

            comm.barrier();

            if (comm.is_root()) {
                moonolith::logger() << "ConvertTransferAlgorithm:apply(...)/Vol2Surf end" << std::endl;
            }

            return true;
        }
    };

    template <>
    class ApplyAux<1, 1, 0> final {
    public:
        using MeshT = moonolith::Mesh<double, 1>;
        using FunctionSpaceT = moonolith::FunctionSpace<MeshT>;

        static bool apply(const FunctionSpaceT &master,
                          const FunctionSpaceT &slave,
                          const TransferOptions &opts,
                          TransferDataT &data) {
            assert(false && "invalid dimensions");
            return false;
        }
    };

    template <>
    class ApplyAux<2, 2, 1> final : public ApplyAuxLowerDim<2> {};
    template <>
    class ApplyAux<3, 3, 2> final : public ApplyAuxLowerDim<3> {};

    template <int Dim>
    class ConvertTransferAlgorithm {
    public:
        using MeshT = moonolith::Mesh<double, Dim>;
        using FunctionSpaceT = moonolith::FunctionSpace<MeshT>;

        template <int DimMaster, int DimSlave>
        static bool apply_aux(const FunctionSpaceT &master,
                              const FunctionSpaceT &slave,
                              const TransferOptions &opts,
                              TransferDataT &data) {
            return ApplyAux<Dim, DimMaster, DimSlave>::apply(master, slave, opts, data);
        }

        static bool apply(const FunctionSpaceT &master,
                          const FunctionSpaceT &slave,
                          const TransferOptions &opts,
                          TransferDataT &data) {
            if (slave.mesh().manifold_dim() < Dim) {
                if (master.mesh().manifold_dim() < Dim) {
                    // surf to surf
                    assert(Dim > 1);
                    return apply_aux<Dim - 1, Dim - 1>(master, slave, opts, data);
                } else {
                    // Vol 2 surf
                    assert(Dim > 1);
                    return apply_aux<Dim, Dim - 1>(master, slave, opts, data);
                }
            } else {
                // vol 2 vol
                apply_aux<Dim, Dim>(master, slave, opts, data);
            }

            return true;
        }

        static bool surface_apply(const libMesh::MeshBase &from_mesh,
                                  const libMesh::DofMap &from_dofs,
                                  const libMesh::MeshBase &to_mesh,
                                  const libMesh::DofMap &to_dofs,
                                  const TransferOptions &opts,
                                  TransferDataT &data) {
            moonolith::Communicator comm = from_mesh.comm().get();
            auto master_mesh = std::make_shared<MeshT>(comm);
            auto slave_mesh = std::make_shared<MeshT>(comm);

            FunctionSpaceT master(master_mesh), slave(slave_mesh);

            extract_trace_space(from_mesh, from_dofs, opts.from_var_num, master);
            extract_trace_space(to_mesh, to_dofs, opts.to_var_num, slave);

            return apply(master, slave, opts, data);
        }

        static bool surface_apply(const libMesh::MeshBase &lm_mesh,
                                  const libMesh::DofMap &lm_dofs,
                                  const TransferOptions &opts,
                                  TransferDataT &data) {
            moonolith::Communicator comm = moonolith::Communicator(lm_mesh.comm().get(), false).duplicate();
            auto mesh = std::make_shared<MeshT>(comm);

            FunctionSpaceT space(mesh);
            extract_trace_space(lm_mesh, lm_dofs, opts.from_var_num, space);

            comm.barrier();

            if (comm.is_root()) {
                moonolith::logger() << "ConvertTransferAlgorithm:surface_apply(...) begin" << std::endl;
            }

            static const int ManifoldDim = moonolith::StaticMax<Dim - 1, 1>::value;

            moonolith::ParL2Transfer<double, Dim, ManifoldDim, ManifoldDim> assembler(comm);

            if (!assembler.assemble(space, opts.tags, 1e-8)) {
                return false;
            }

            prepare_data(opts, assembler, data);

            comm.barrier();

            if (comm.is_root()) {
                moonolith::logger() << "ConvertTransferAlgorithm:surface_apply(...) end" << std::endl;
            }

            return true;
        }

        static bool apply(const libMesh::MeshBase &from_mesh,
                          const libMesh::DofMap &from_dofs,
                          const libMesh::MeshBase &to_mesh,
                          const libMesh::DofMap &to_dofs,
                          const TransferOptions &opts,
                          TransferDataT &data) {
            moonolith::Communicator comm = moonolith::Communicator(from_mesh.comm().get(), false).duplicate();
            auto master_mesh = std::make_shared<MeshT>(comm);
            auto slave_mesh = std::make_shared<MeshT>(comm);

            FunctionSpaceT master(master_mesh), slave(slave_mesh);

            convert_libmesh_to_moonolith(from_mesh, from_dofs, opts.from_var_num, master);
            convert_libmesh_to_moonolith(to_mesh, to_dofs, opts.to_var_num, slave);

            return apply(master, slave, opts, data);
        }

        template <class Transfer>
        static void prepare_data_with_covering_check(const TransferOptions &opts,
                                                     const FunctionSpaceT &in_slave,
                                                     const FunctionSpaceT &slave,
                                                     Transfer &t,
                                                     TransferDataT &data) {
            {
                USparseMatrix permutation;

                if (opts.n_var > 1) {
                    make_tensorize_permutation(opts.n_var, slave, in_slave, permutation);
                } else {
                    make_permutation(slave, in_slave, permutation);
                }

                TransferDataT data_ew;

                auto &B_ew = *data_ew.B;
                auto &D_ew = *data_ew.D;
                auto &Q_ew = *data_ew.Q;
                // auto &T_ew = *data_ew.T;

                convert_matrix(t.buffers.B.get(), B_ew);
                convert_matrix(t.buffers.D.get(), D_ew);
                convert_matrix(t.buffers.Q.get(), Q_ew);

                data_ew.permute(permutation, data);
            }

            auto &B = *data.B;
            auto &D = *data.D;
            auto &Q = *data.Q;
            auto &T = *data.T;

            if (!empty(Q)) {
                m_utopia_warning_once("using sum(D, 1) instead of diag(D)");
                UVector d_inv = sum(D, 1);

                e_pseudo_inv(d_inv, d_inv, 1e-12);

                USparseMatrix D_tilde_inv = diag(d_inv);
                USparseMatrix T_temp = D_tilde_inv * B;

                if (opts.n_var == 1) {
                    T = Q * T_temp;
                } else {
                    USparseMatrix T_x = Q * T_temp;
                    tensorize(T_x, opts.n_var, T);
                }
            }

            // write("T.m", T);
        }

        static bool apply_with_covering_check(const FunctionSpaceT &master,
                                              const FunctionSpaceT &in_slave,
                                              const TransferOptions &opts,
                                              TransferDataT &data) {
            FunctionSpaceT slave;
            in_slave.separate_dofs(slave);

            moonolith::Communicator comm = master.mesh().comm();

            comm.barrier();

            if (comm.is_root()) {
                moonolith::logger() << "ConvertTransferAlgorithm:apply_with_covering_check(...) begin" << std::endl;
            }

            if (slave.mesh().manifold_dim() > Dim) {
                assert(Dim > 1);
                moonolith::ParL2Transfer<double, Dim, Dim, moonolith::StaticMax<Dim - 1, 1>::value> assembler(comm);
                assembler.remove_incomplete_intersections(true);

                if (opts.tags.empty()) {
                    // change assemble with covering check
                    if (!assembler.assemble(master, slave)) {
                        return false;
                    }
                } else {
                    // change assemble with covering check
                    if (!assembler.assemble(master, slave, opts.tags)) {
                        return false;
                    }
                }

                prepare_data_with_covering_check(opts, in_slave, slave, assembler, data);
            } else {
                moonolith::ParL2Transfer<double, Dim, Dim, Dim> assembler(comm);
                assembler.remove_incomplete_intersections(true);

                if (opts.tags.empty()) {
                    // change assemble with covering check
                    if (!assembler.assemble(master, slave)) {
                        return false;
                    }
                } else {
                    // change assemble with covering check
                    if (!assembler.assemble(master, slave, opts.tags)) {
                        return false;
                    }
                }

                prepare_data_with_covering_check(opts, in_slave, slave, assembler, data);
            }

            comm.barrier();
            if (comm.is_root()) {
                moonolith::logger() << "ConvertTransferAlgorithm:apply_with_covering_check(...) end" << std::endl;
            }

            return true;
        }

        static bool apply_with_covering_check(const libMesh::MeshBase &from_mesh,
                                              const libMesh::DofMap &from_dofs,
                                              const libMesh::MeshBase &to_mesh,
                                              const libMesh::DofMap &to_dofs,
                                              const TransferOptions &opts,
                                              TransferDataT &data) {
            moonolith::Communicator comm = moonolith::Communicator(from_mesh.comm().get(), false).duplicate();

            auto master_mesh = std::make_shared<MeshT>(comm);
            auto slave_mesh = std::make_shared<MeshT>(comm);

            FunctionSpaceT master(master_mesh), slave(slave_mesh);

            convert_libmesh_to_moonolith(from_mesh, from_dofs, opts.from_var_num, master);
            convert_libmesh_to_moonolith(to_mesh, to_dofs, opts.to_var_num, slave);

            return apply_with_covering_check(master, slave, opts, data);
        }
    };

    bool NewTransferAssembler::surface_assemble(const std::shared_ptr<MeshBase> &from_mesh,
                                                const std::shared_ptr<DofMap> &from_dofs,
                                                const std::shared_ptr<MeshBase> &to_mesh,
                                                const std::shared_ptr<DofMap> &to_dofs,
                                                const TransferOptions &opts) {
        auto spatial_dim = to_mesh->spatial_dimension();
        bool has_intersection = false;

        // std::cout << "handle_adaptive_refinement_ : " << handle_adaptive_refinement_ << std::endl;

        if (handle_adaptive_refinement_) {
            Adaptivity a;
            a.constraint_matrix(
                *from_mesh, *from_dofs, *data.constraint_matrix_from, *data.post_constraint_matrix_from);
            a.constraint_matrix(*to_mesh, *to_dofs, *data.constraint_matrix_to, *data.post_constraint_matrix_to);

            // disp("from");
            // disp(*data.constraint_matrix_from);
            // disp("to");
            // disp(*data.constraint_matrix_to);
        }

        if (spatial_dim == 1) {
            // if(remove_incomplete_intersections_) {
            //     has_intersection = ConvertTransferAlgorithm<1>::apply_with_covering_check(*from_mesh, *from_dofs,
            //     *to_mesh, *to_dofs, opts, data);
            // } else {
            has_intersection =
                ConvertTransferAlgorithm<1>::surface_apply(*from_mesh, *from_dofs, *to_mesh, *to_dofs, opts, data);
            // }
        } else if (spatial_dim == 2) {
            // if(remove_incomplete_intersections_) {
            //     has_intersection = ConvertTransferAlgorithm<2>::apply_with_covering_check(*from_mesh, *from_dofs,
            //     *to_mesh, *to_dofs, opts, data);
            // } else {
            has_intersection =
                ConvertTransferAlgorithm<2>::surface_apply(*from_mesh, *from_dofs, *to_mesh, *to_dofs, opts, data);
            // }
        } else if (spatial_dim == 3) {
            // if(remove_incomplete_intersections_) {
            //     has_intersection = ConvertTransferAlgorithm<3>::apply_with_covering_check(*from_mesh, *from_dofs,
            //     *to_mesh, *to_dofs, opts, data);
            // } else {
            has_intersection =
                ConvertTransferAlgorithm<3>::surface_apply(*from_mesh, *from_dofs, *to_mesh, *to_dofs, opts, data);
            // }
        }

        // write("T.m", *data.T);

        return has_intersection;
    }

    bool NewTransferAssembler::assemble(const std::shared_ptr<MeshBase> &from_mesh,
                                        const std::shared_ptr<DofMap> &from_dofs,
                                        const std::shared_ptr<MeshBase> &to_mesh,
                                        const std::shared_ptr<DofMap> &to_dofs,
                                        const TransferOptions &opts) {
        // std::cout << "handle_adaptive_refinement_ : " << handle_adaptive_refinement_ << std::endl;

        if (handle_adaptive_refinement_) {
            Adaptivity a;
            a.constraint_matrix(
                *from_mesh, *from_dofs, *data.constraint_matrix_from, *data.post_constraint_matrix_from);
            a.constraint_matrix(*to_mesh, *to_dofs, *data.constraint_matrix_to, *data.post_constraint_matrix_to);

            // disp("from");
            // disp(*data.constraint_matrix_from);
            // disp("to");
            // disp(*data.constraint_matrix_to);
        }

        if (use_convert_transfer_) {
            auto spatial_dim = to_mesh->spatial_dimension();
            bool has_intersection = false;

            if (spatial_dim == 1) {
                if (remove_incomplete_intersections_) {
                    has_intersection = ConvertTransferAlgorithm<1>::apply_with_covering_check(
                        *from_mesh, *from_dofs, *to_mesh, *to_dofs, opts, data);
                } else {
                    has_intersection =
                        ConvertTransferAlgorithm<1>::apply(*from_mesh, *from_dofs, *to_mesh, *to_dofs, opts, data);
                }
            } else if (spatial_dim == 2) {
                if (remove_incomplete_intersections_) {
                    has_intersection = ConvertTransferAlgorithm<2>::apply_with_covering_check(
                        *from_mesh, *from_dofs, *to_mesh, *to_dofs, opts, data);
                } else {
                    has_intersection =
                        ConvertTransferAlgorithm<2>::apply(*from_mesh, *from_dofs, *to_mesh, *to_dofs, opts, data);
                }
            } else if (spatial_dim == 3) {
                if (remove_incomplete_intersections_) {
                    has_intersection = ConvertTransferAlgorithm<3>::apply_with_covering_check(
                        *from_mesh, *from_dofs, *to_mesh, *to_dofs, opts, data);
                } else {
                    has_intersection =
                        ConvertTransferAlgorithm<3>::apply(*from_mesh, *from_dofs, *to_mesh, *to_dofs, opts, data);
                }
            }

            // write("T.m", *data.T);

            return has_intersection;

        } else {
            LibMeshFunctionSpaceAdapter from_adapter, to_adapter;
            from_adapter.init(from_mesh, *from_dofs, opts.from_var_num);
            to_adapter.init(to_mesh, *to_dofs, opts.to_var_num);

            auto spatial_dim = to_mesh->spatial_dimension();

            bool has_intersection = false;

            if (spatial_dim == 1) {
                has_intersection = TransferAlgorithm<1>::apply(from_adapter, to_adapter, opts, data);
            } else if (spatial_dim == 2) {
                has_intersection = TransferAlgorithm<2>::apply(from_adapter, to_adapter, opts, data);
            } else if (spatial_dim == 3) {
                has_intersection = TransferAlgorithm<3>::apply(from_adapter, to_adapter, opts, data);
            }

            return has_intersection;
        }
    }

    bool NewTransferAssembler::assemble(const MeshBase &from_mesh,
                                        const DofMap &from_dofs,
                                        const MeshBase &to_mesh,
                                        const DofMap &to_dofs,
                                        const TransferOptions &opts) {
        auto spatial_dim = to_mesh.spatial_dimension();
        bool has_intersection = false;

        // std::cout << "handle_adaptive_refinement_ : " << handle_adaptive_refinement_ << std::endl;

        if (handle_adaptive_refinement_) {
            Adaptivity a;
            // a.constraint_matrix(*from_mesh, *from_dofs, opts.from_var_num, *data.constraint_matrix_from,
            // *data.constraint_matrix_from_2); a.constraint_matrix(*to_mesh,   *to_dofs,   opts.to_var_num,
            // *data.constraint_matrix_to,*data.constraint_matrix_to_2);
            a.constraint_matrix(from_mesh, from_dofs, *data.constraint_matrix_from, *data.post_constraint_matrix_from);
            a.constraint_matrix(to_mesh, to_dofs, *data.constraint_matrix_to, *data.post_constraint_matrix_to);

            // write("from.m", *data.constraint_matrix_from);
            // write("to.m", *data.constraint_matrix_to);

            // write("from_post.m", *data.post_constraint_matrix_from);
            // write("to_post.m", *data.post_constraint_matrix_to);

            // write("from_mesh.m", from_mesh);
            // write("to_mesh.m", *to_mesh);
        }

        if (spatial_dim == 1) {
            if (remove_incomplete_intersections_) {
                has_intersection = ConvertTransferAlgorithm<1>::apply_with_covering_check(
                    from_mesh, from_dofs, to_mesh, to_dofs, opts, data);
            } else {
                has_intersection =
                    ConvertTransferAlgorithm<1>::apply(from_mesh, from_dofs, to_mesh, to_dofs, opts, data);
            }
        } else if (spatial_dim == 2) {
            if (remove_incomplete_intersections_) {
                has_intersection = ConvertTransferAlgorithm<2>::apply_with_covering_check(
                    from_mesh, from_dofs, to_mesh, to_dofs, opts, data);
            } else {
                has_intersection =
                    ConvertTransferAlgorithm<2>::apply(from_mesh, from_dofs, to_mesh, to_dofs, opts, data);
            }
        } else if (spatial_dim == 3) {
            if (remove_incomplete_intersections_) {
                has_intersection = ConvertTransferAlgorithm<3>::apply_with_covering_check(
                    from_mesh, from_dofs, to_mesh, to_dofs, opts, data);
            } else {
                has_intersection =
                    ConvertTransferAlgorithm<3>::apply(from_mesh, from_dofs, to_mesh, to_dofs, opts, data);
            }
        }

        return has_intersection;
    }

    bool NewTransferAssembler::surface_assemble(const MeshBase &mesh, const DofMap &dofs, const TransferOptions &opts) {
        if (handle_adaptive_refinement_) {
            Adaptivity a;
            a.constraint_matrix(mesh, dofs, *data.constraint_matrix_from, *data.post_constraint_matrix_from);
            a.constraint_matrix(mesh, dofs, *data.constraint_matrix_to, *data.post_constraint_matrix_to);

            // disp("from");
            // disp(*data.constraint_matrix_from);
            // disp("to");
            // disp(*data.constraint_matrix_to);
        }

        auto spatial_dim = mesh.spatial_dimension();
        bool has_intersection = false;

        if (spatial_dim == 1) {
            has_intersection = ConvertTransferAlgorithm<1>::surface_apply(mesh, dofs, opts, data);
        } else if (spatial_dim == 2) {
            has_intersection = ConvertTransferAlgorithm<2>::surface_apply(mesh, dofs, opts, data);
        } else if (spatial_dim == 3) {
            has_intersection = ConvertTransferAlgorithm<3>::surface_apply(mesh, dofs, opts, data);
        }

        return has_intersection;
    }
}  // namespace utopia

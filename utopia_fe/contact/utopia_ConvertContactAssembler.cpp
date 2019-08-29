#include "utopia_ConvertContactAssembler.hpp"

#include "utopia_DualBasis.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_MatrixInserter.hpp"
#include "utopia_ZeroRowsToIdentity.hpp"
#include "utopia_NormalizeRows.hpp"
#include "utopia_LibMeshShape.hpp"
#include "utopia_SurfaceQuadratureConverter.hpp"
#include "utopia_OrthogonalTransformation.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_TransferUtils.hpp"
#include "utopia_moonolith_permutations.hpp"

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
#include "moonolith_par_contact.hpp"

#include <vector>
#include <memory>

#include "libmesh/mesh_refinement.h"

namespace utopia {

    inline static bool check_op(const USparseMatrix &T)
    {
        UVector sum_T = sum(T, 1);
        
        bool ok = true;
        each_read(sum_T, [&ok](const SizeType i, const double val) {
            if(!approxeq(val, 0.0, 1e-1) && !approxeq(val, 1.0, 1e-1)) {
                std::cerr << i << "] " << val << "\n";
                ok = false;
            }
        });

        if(!ok) {
            write("T_bad.m", T);
        }
        
        return ok;
    }


    bool ConvertContactTensors::write()
    {
        if(!empty(B)) 
        {
            B.implementation().set_name("b");

            utopia::write("B.m", B);
        }

        if(!empty(D)) 
        {
            D.implementation().set_name("d");

            utopia::write("D.m", D);
        }

        if(!empty(Q)) 
        {
            Q.implementation().set_name("q");

            utopia::write("Q.m", Q);
        }

        if(!empty(Q_inv)) 
        {
            Q_inv.implementation().set_name("q_inv");

            utopia::write("Q_inv.m", Q_inv);
        }

        if(!empty(T)) 
        {
            T.implementation().set_name("t");

            utopia::write("T.m", T);
        }

        if(!empty(inv_mass_vector)) 
        {
            inv_mass_vector.implementation().set_name("imv");

            utopia::write("inv_mass_vector.m", inv_mass_vector);
        }

        if(!empty(weighted_gap)) 
        {
            weighted_gap.implementation().set_name("wg");

            utopia::write("weighted_gap.m", weighted_gap);
        }

        if(!empty(gap)) 
        {
            gap.implementation().set_name("g");

            utopia::write("gap.m", gap);
        }

        if(!empty(weighted_normal)) 
        {
            weighted_normal.implementation().set_name("wn");

            utopia::write("weighted_normal.m", weighted_normal);
        }

        if(!empty(normal)) 
        {
            normal.implementation().set_name("n");

            utopia::write("normal.m", normal);
        }


        if(!empty(is_contact)) 
        {
            is_contact.implementation().set_name("ic");

            utopia::write("is_contact.m", is_contact);
        }


        if(!empty(is_glue)) 
        {
            is_glue.implementation().set_name("ig");

            utopia::write("is_glue.m", is_glue);
        }

// orthogonal_trafo.implementation().set_name("o");
// utopia::write("O.m", orthogonal_trafo);

// complete_transformation.implementation().set_name("ct");
// utopia::write("complete_transformation.m", complete_transformation);
        return true;
    }
    


    class ConvertContactBuffers {
    public:

        template<int Dim>
        using SpaceT = moonolith::FunctionSpace<moonolith::Mesh<double, Dim>>;

        template<int Dim>
        void finalize(
            const moonolith::ContactBuffers<double> &buffers,
            const SpaceT<Dim> &element_wise_space,
            const SpaceT<Dim> &node_wise_space)
        {
            convert_matrix(buffers.B,          element_wise.B);
            convert_matrix(buffers.D,          element_wise.D);
            convert_matrix(buffers.Q,          element_wise.Q);
            convert_matrix(buffers.Q_inv,      element_wise.Q_inv);
            
            convert_tensor(buffers.gap,        element_wise.weighted_gap);
            convert_tensor(buffers.normal,     element_wise.weighted_normal);
            convert_tensor(buffers.is_glue,    element_wise.is_glue);
            convert_tensor(buffers.is_contact, element_wise.is_contact);

            convert_to_node_wise(element_wise_space, element_wise, node_wise_space, *node_wise);
        }

        template<int Dim>
        static void normalize_and_build_orthgonal_trafo(const SpaceT<Dim> &node_wise_space, ConvertContactTensors &node_wise)
        {
            node_wise.orthogonal_trafo = local_sparse(node_wise_space.dof_map().n_local_dofs(), node_wise_space.dof_map().n_local_dofs(), Dim);

            moonolith::HouseholderTransformation<double, Dim> H;
            auto &normal = node_wise.normal;

            auto r = range(normal);

            ReadAndWrite<UVector> rw_normal(normal);
            Read<UVector> r_is_c(node_wise.is_contact);
            Write<USparseMatrix> w_ot(node_wise.orthogonal_trafo, utopia::LOCAL);

            moonolith::Vector<double, Dim> n;
            for(auto i = r.begin(); i < r.end(); i += Dim) {
                const bool is_contact = node_wise.is_contact.get(i);

                if(is_contact) {
                    
                    for(int d = 0; d < Dim; ++d) {
                        n[d] = normal.get(i + d);
                    }

                    auto len = length(n);
                    assert(len > 0.0);

                    n /= len;

                    for(int d = 0; d < Dim; ++d) {
                        normal.set(i + d, n[d]);
                    }

                    n.x -= 1.0;

                    len = length(n);

                    if(approxeq(len, 0.0, 1e-15)) {
                        H.identity();
                    } else {
                        n /= len;
                        H.init(n);
                    }

                    assert(approxeq(std::abs(measure(H)), 1.0, 1e-8));

                    for(int d1 = 0; d1 < Dim; ++d1) {
                        for(int d2 = 0; d2 < Dim; ++d2) {
                            node_wise.orthogonal_trafo.set(i + d1, i + d2, H(d1, d2));
                        }
                    }

                } else {

                    for(int d1 = 0; d1 < Dim; ++d1) {
                        node_wise.orthogonal_trafo.set(i + d1, i + d1, 1.0);
                    }

                }
            }

        }

        template<int Dim>
        static void convert_to_node_wise(
            const SpaceT<Dim> &elem_wise_space, ConvertContactTensors &elem_wise,
            const SpaceT<Dim> &node_wise_space, ConvertContactTensors &node_wise)
        {
            USparseMatrix perm, vector_perm;

            make_permutation(
                elem_wise_space,
                node_wise_space,
                perm
            );

            make_vector_permutation(
                Dim,
                elem_wise_space,
                node_wise_space,
                vector_perm
            );

            node_wise.weighted_gap = perm * elem_wise.weighted_gap;
            node_wise.is_glue      = perm * elem_wise.is_glue;
            node_wise.is_contact   = perm * elem_wise.is_contact;

            node_wise.weighted_normal = vector_perm * elem_wise.weighted_normal;

            USparseMatrix B_x = perm * elem_wise.B * transpose(perm);
            USparseMatrix D_x = perm * elem_wise.D * transpose(perm);
            USparseMatrix Q_x = perm * elem_wise.Q * transpose(perm);

            USparseMatrix Q_inv_x = perm * elem_wise.Q_inv * transpose(perm);

            normalize_rows(Q_x,     1e-15);
            normalize_rows(Q_inv_x, 1e-15);

            {
                each_transform(node_wise.is_contact, node_wise.is_contact, [&node_wise](const SizeType i, const double value) -> double
                {
                    if(value > 0.0) {
                        return 1.0;
                    } else {
                        return 0.0;
                    }
                });

                each_transform(node_wise.is_glue, node_wise.is_glue, [&node_wise](const SizeType i, const double value) -> double
                {
                    if(value > 0.0) {
                        return 1.0;
                    } else {
                        return 0.0;
                    }
                });
            }

            node_wise.inv_mass_vector = sum(D_x, 1);

            e_pseudo_inv(node_wise.inv_mass_vector, node_wise.inv_mass_vector, 1e-15);
            USparseMatrix D_inv_x = diag(node_wise.inv_mass_vector);

            USparseMatrix T_temp_x = D_inv_x * B_x;
            USparseMatrix T_x      = Q_x * T_temp_x;

            tensorize(B_x,     Dim, node_wise.B);
            tensorize(D_x,     Dim, node_wise.D);
            tensorize(T_x,     Dim, node_wise.T);
            tensorize(Q_x,     Dim, node_wise.Q);
            tensorize(Q_inv_x, Dim, node_wise.Q_inv);
            
            tensorize(Dim, node_wise.inv_mass_vector);
            tensorize(Dim, node_wise.is_glue);
            
            normalize_rows(node_wise.T, 1e-15);

            assert(check_op(node_wise.T));

            node_wise.T += local_identity(local_size(node_wise.T));
            
            node_wise.gap    = node_wise.Q * e_mul(node_wise.inv_mass_vector, node_wise.weighted_gap);
            node_wise.normal = node_wise.Q * e_mul(node_wise.inv_mass_vector, node_wise.weighted_normal);

            normalize_and_build_orthgonal_trafo(node_wise_space, node_wise);
            node_wise.complete_transformation = node_wise.T * node_wise.orthogonal_trafo;

            {
                static const double LARGE_VALUE = 1e6;
                Write<UVector> w(node_wise.gap);
                each_read(node_wise.is_contact, [&node_wise](const SizeType i, const double value)
                {
                    if(value == 0.0) {
                        node_wise.gap.set(i, LARGE_VALUE);
                    }
                });
            }

            // node_wise.write();
        }

        ConvertContactTensors element_wise;
        std::shared_ptr<ConvertContactTensors> node_wise;

        ConvertContactBuffers(MPI_Comm comm) {
            node_wise = std::make_shared<ConvertContactTensors>();
        }

    private:
        ConvertContactBuffers(const ConvertContactBuffers &) {}
    };

    template<int Dim>
    class ConvertContactAlgorithm {
    public:
        using MeshT  = moonolith::Mesh<double, Dim>;
        using SpaceT = moonolith::FunctionSpace<MeshT>;

        static bool apply(const ContactParams &params,
                          const std::shared_ptr<ElementBlackList> &black_list,
                          libMesh::MeshBase &mesh,
                          libMesh::DofMap &dof_map,
                          ConvertContactBuffers &contact_data)
        {
            using AlogrithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, LibMeshFunctionSpaceAdapter>;
            using Adapter    = typename AlogrithmT::Adapter;

            moonolith::Communicator comm = mesh.comm().get();

            auto mesh_ptr = std::make_shared<MeshT>(comm);
            SpaceT space(mesh_ptr);

            extract_trace_space(mesh, dof_map, params.variable_number, space);
            // moonolith::MatlabScripter script;
            // mesh_ptr->draw(script);
            // script.save("contact.m");

            SpaceT elem_wise_space;
            space.separate_dofs(elem_wise_space);

            moonolith::ParContact<double, Dim> par_contact(comm, Dim == 2);

            if(par_contact.assemble(
                params.contact_pair_tags,
                elem_wise_space,
                params.side_set_search_radius,
                params.is_glue)
            ) {
                contact_data.finalize(par_contact.buffers, elem_wise_space, space);
                return true;
            } else {
                return false;
            }
        }


        //without dof separation
        // static bool apply(const ContactParams &params,
        //                   const std::shared_ptr<ElementBlackList> &black_list,
        //                   libMesh::MeshBase &mesh,
        //                   libMesh::DofMap &dof_map,
        //                   ConvertContactBuffers &contact_data)
        // {
        //     using AlogrithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, LibMeshFunctionSpaceAdapter>;
        //     using Adapter    = typename AlogrithmT::Adapter;

        //     moonolith::Communicator comm = mesh.comm().get();

        //     auto mesh_ptr = std::make_shared<MeshT>(comm);
        //     SpaceT space(mesh_ptr);

        //     extract_trace_space(mesh, dof_map, params.variable_number, space);
        //     // moonolith::MatlabScripter script;
        //     // mesh_ptr->draw(script);
        //     // script.save("contact.m");

        //     // SpaceT elem_wise_space;
        //     // space.separate_dofs(elem_wise_space);

        //     moonolith::ParContact<double, Dim> par_contact(comm, Dim == 2);

        //     if(par_contact.assemble(
        //         params.contact_pair_tags,
        //         space,
        //         params.side_set_search_radius,
        //         params.is_glue)
        //     ) {
        //         contact_data.finalize(par_contact.buffers, space, space);
        //         return true;
        //     } else {
        //         return false;
        //     }
        // }
    };

    bool ConvertContactAssembler::assemble(
                                    libMesh::MeshBase &mesh,
                                    libMesh::DofMap &dof_map,
                                    const ContactParams &params)
    {
        Chrono overall_time;
        overall_time.start();

        ConvertContactBuffers contact_data(mesh.comm().get());

        const int spatial_dim = mesh.spatial_dimension();

        has_contact_ = false;
        if(spatial_dim == 2) {
            has_contact_ = ConvertContactAlgorithm<2>::apply(params, black_list_, mesh, dof_map, contact_data);
        } else if(spatial_dim == 3) {
            has_contact_ = ConvertContactAlgorithm<3>::apply(params, black_list_, mesh, dof_map, contact_data);
        }

        if(has_contact_) {
            contact_tensors_ = contact_data.node_wise;
        } else {
            //init default
            init_no_contact(
                mesh,
                dof_map);
        }

        has_glue_ = has_contact_ && !params.is_glue->empty();

        if(has_glue_) {
            double n_glued = sum(contact_tensors_->is_glue);

            if(n_glued == 0.) {
                has_glue_ = false;
            }
        }

        overall_time.stop();
        std::cout << "ConvertContactAssembler::assemble: " << overall_time << std::endl;
        return has_contact_;
    }

    void ConvertContactAssembler::couple(const UVector &in, UVector &out) const
    {
        assert(contact_tensors_);
        out = transpose(contact_tensors_->complete_transformation) * in;
    }

    void ConvertContactAssembler::uncouple(const UVector &in, UVector &out) const
    {
        assert(contact_tensors_);
        out = contact_tensors_->complete_transformation * in;
    }

    void ConvertContactAssembler::couple(const USparseMatrix &in, USparseMatrix &out) const
    {
        assert(contact_tensors_);

        const auto &T = contact_tensors_->complete_transformation;

        out = transpose(T) * in * T;
    }

    const UVector &ConvertContactAssembler::gap() const
    {
        assert(contact_tensors_);
        return contact_tensors_->gap;
    }

    UVector &ConvertContactAssembler::gap()
    {
        assert(contact_tensors_);
        return contact_tensors_->gap;
    }

    bool ConvertContactAssembler::init_no_contact(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map)
    {

        if(!contact_tensors_) {
            contact_tensors_ = std::make_shared<ConvertContactTensors>();
        }

        auto n_local_dofs = dof_map.n_local_dofs();

        contact_tensors_->gap = local_values(n_local_dofs, 100000000);
        contact_tensors_->weighted_gap = local_values(n_local_dofs, 100000000);
        contact_tensors_->normal = local_zeros(n_local_dofs);

        contact_tensors_->inv_mass_vector = local_values(n_local_dofs, 1.);
        contact_tensors_->is_contact = local_zeros(n_local_dofs);

        contact_tensors_->T = local_identity(n_local_dofs, n_local_dofs);
        contact_tensors_->orthogonal_trafo  = local_identity(n_local_dofs, n_local_dofs);
        contact_tensors_->complete_transformation = local_identity(n_local_dofs, n_local_dofs);
        has_contact_ = false;
        has_glue_ = false;
        return true;
    }

    void ConvertContactAssembler::remove_mass(const UVector &in, UVector &out) const
    {
        if(!empty(contact_tensors_->Q)) {
            out = contact_tensors_->Q * e_mul(contact_tensors_->inv_mass_vector, in);
            return;
        }

        out = e_mul(contact_tensors_->inv_mass_vector, in);
    }

    void ConvertContactAssembler::read(Input &in)
    {
        in.get("black-list", [this](Input &in) {
            black_list_ = std::make_shared<ElementBlackList>(true);
            black_list_->read(in);
        });
    }
}


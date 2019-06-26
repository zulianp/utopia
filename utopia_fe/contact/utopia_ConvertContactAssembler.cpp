// #include "utopia_ConvertContactAssembler.hpp"

// #include "utopia_DualBasis.hpp"
// #include "utopia_libmesh_NonLinearFEFunction.hpp"
// #include "utopia_LibMeshBackend.hpp"
// #include "utopia_MatrixInserter.hpp"
// #include "utopia_ZeroRowsToIdentity.hpp"
// #include "utopia_NormalizeRows.hpp"
// #include "utopia_LibMeshShape.hpp"
// #include "utopia_SurfaceQuadratureConverter.hpp"
// #include "utopia_OrthogonalTransformation.hpp"
// #include "utopia_ElementWisePseudoInverse.hpp"

// #include "moonolith_affine_transform.hpp"
// #include "moonolith_contact.hpp"
// #include "moonolith_sparse_matrix.hpp"
// #include "moonolith_redistribute.hpp"
// #include "moonolith_assign_functions.hpp"
// #include "moonolith_elem_shape.hpp"
// #include "moonolith_elem_triangle.hpp"
// #include "moonolith_elem_quad.hpp"
// #include "moonolith_elem_segment.hpp"
// #include "moonolith_matlab_scripter.hpp"
// #include "moonolith_par_contact.hpp"

// #include <vector>
// #include <memory>

// #include "libmesh/mesh_refinement.h"

// namespace utopia {

//     class ConvertContactBuffers {
//     public:

//         inline MPI_Comm comm() const
//         {
//             return B.comm.get();
//         }


//         template<int Dim>
//         void finalize(
//             const moonolith::FunctionSpace<moonolith::Mesh<double, Dim>> &element_wise_space,
//             const moonolith::FunctionSpace<moonolith::Mesh<double, Dim>> &dof_wise_space)
//         {
//             const SizeType n_local_elems = element_wise_space.n_elements();
//             const SizeType n_local_dofs  = element_wise_space.n_local_dofs();
//             const SizeType spatial_dim   = Dim;


//             area.finalize(n_local_elems);
//             area.fill(element_wise.area);

//             auto r = range(volumes);
//             assert(r.extent() == n_local_elems);

//             std::vector<bool> remove(adapter.n_local_dofs(), false);
//             auto cr = adapter.permutation()->implementation().col_range();

//             element_wise.is_contact = local_zeros(adapter.n_local_dofs());

//             UVector elem_to_transform;

//             const bool must_assemble_trafo = use_biorth && adapter.fe_type(0).order == 2;
//             elem_to_transform = local_zeros(r.extent());

//             {
//                 Read<UVector> rv(volumes), ra(element_wise.area);
//                 Write<UVector> wic(element_wise.is_contact), wett(elem_to_transform);

//                 for(auto i = r.begin(); i < r.end(); ++i) {
//                     const auto a = element_wise.area.get(i);

//                     if(a > 0.0) {
//                         auto ratio = a/volumes.get(i);

//                         if(!approxeq(ratio, 1.0, 1e-2)) {
//                             remove[i - r.begin()] = true;

//                             const auto &dofs = adapter.element_dof_map()[i - r.begin()].dofs;

//                             for(auto d : dofs) {
//                                 remove[d - cr.begin()] = true;
//                             }

//                             // std::cout << "=====================================\n";
//                             // std::cout << i << ") removed. " << int(ratio * 100) << "% of slave volume" << std::endl;

//                         } else {
//                             // std::cout << "=====================================\n";
//                             // std::cout << i << ") " << volumes.get(i) << " == " << element_wise.area.get(i) << std::endl;

//                             const auto &dofs = adapter.element_dof_map()[i - r.begin()].dofs;

//                             for(auto d : dofs) {
//                                 element_wise.is_contact.set(d, 1.0);
//                             }

//                             if(must_assemble_trafo) {
//                                 elem_to_transform.set(i, 1.);
//                             }
//                         }
//                     }
//                 }
//             }

//             B.finalize(n_local_dofs, n_local_dofs);
//             D.finalize(n_local_dofs, n_local_dofs);
//             gap.finalize(n_local_dofs);
//             normal.finalize(n_local_dofs * spatial_dim);
//             is_glue.finalize(n_local_dofs);

//             B.fill(remove, element_wise.B_x);
//             D.fill(remove, element_wise.D_x);
//             gap.fill(remove, element_wise.weighted_gap);
//             normal.fill(element_wise.weighted_normal);
//             is_glue.fill(element_wise.is_glue);

//             tensor_prod_with_identity(element_wise.B_x, spatial_dim, element_wise.B);
//             tensor_prod_with_identity(element_wise.D_x, spatial_dim, element_wise.D);

//             if(must_assemble_trafo) {
//                 DualBasis::build_global_trafo(
//                                               adapter.mesh(),
//                                               adapter.n_local_dofs(),
//                                               adapter.element_dof_map(),
//                                               elem_to_transform,
//                                               alpha,
//                                               element_wise.Q_x,
//                                               false
//                                               );


//                 tensor_prod_with_identity(element_wise.Q_x, spatial_dim, element_wise.Q);


//                 DualBasis::build_global_trafo(
//                                               adapter.mesh(),
//                                               adapter.n_local_dofs(),
//                                               adapter.element_dof_map(),
//                                               elem_to_transform,
//                                               alpha,
//                                               element_wise.Q_x,
//                                               true);


//                 tensor_prod_with_identity(element_wise.Q_x, spatial_dim, element_wise.Q_inv);
//             }

//             double sum_B_x = sum(element_wise.B_x);
//             double sum_D_x = sum(element_wise.D_x);
//             double sum_normal_xyz = sum(element_wise.weighted_normal);

//             assert(adapter.permutation());

//             element_wise.convert(*adapter.permutation(),
//                                  *adapter.vector_permutation(),
//                                  *dof_wise
//                                  );

//             dof_wise->finalize(spatial_dim);


//             static const double LARGE_VALUE = 1e6;

//             Read<UVector> ric(dof_wise->is_contact);
//             each_transform(dof_wise->gap, dof_wise->gap, [&](const SizeType i, const double value) -> double {
//                 const SizeType offset = i;

//                 if(dof_wise->is_contact.get(i) > 0) {
//                     return value;
//                 } else {
//                     return LARGE_VALUE;
//                 }
//             });
//         }

//         MatrixInserter B, D;
//         MatrixInserter gap, normal;
//         MatrixInserter area;
//         MatrixInserter is_glue;
//         ContactTensors element_wise;
//         std::shared_ptr<ContactTensors> dof_wise;

//         ConvertContactBuffers(MPI_Comm comm) : B(comm), D(comm), gap(comm), normal(comm), area(comm), is_glue(comm, false) {
//             dof_wise = std::make_shared<ContactTensors>();
//         }
//     private:
//         ConvertContactBuffers(const ConvertContactBuffers &other)
//         : B(other.comm()), D(other.comm()), gap(other.comm()), normal(other.comm()), area(other.comm()), is_glue(other.comm()) {}
//     };

//     template<int Dim>
//     class ConvertContactAlgorithm {
//     public:
//         using MeshT  = moonolith::Mesh<double, Dim>;
//         using SpaceT = moonolith::FunctionSpace<MeshT>;

//         static bool apply(const ContactParams &params,
//                           const std::shared_ptr<ElementBlackList> &black_list,
//                           libMesh::MeshBase &mesh,
//                           libMesh::DofMap &dof_map,
//                           ConvertContactBuffers &contact_data)
//         {
//             using AlogrithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, LibMeshFunctionSpaceAdapter>;
//             using Adapter    = typename AlogrithmT::Adapter;

//             moonolith::Communicator comm = mesh.comm().get();

//             auto mesh_ptr = std::make_shared<MeshT>(comm);
//             SpaceT space(mesh_ptr);

//             extract_trace_space(mesh, dof_map, params.variable_number, space);

//             moonolith::ParContact<double, Dim> par_contact(comm);

//             if(par_contact.assemble(
//                 params.contact_pair_tags,
//                 space,
//                 params.side_set_search_radius,
//                 params.is_glue)
//             ) {

//             } else {
//                 return false;
//             }
//         }
//     };

//     bool ConvertContactAssembler::assemble(
//                                     libMesh::MeshBase &mesh,
//                                     libMesh::DofMap &dof_map,
//                                     const ContactParams &params)
//     {
//         Chrono overall_time;
//         overall_time.start();

//         ConvertContactBuffers contact_data(mesh.comm().get());

//         const int spatial_dim = mesh.spatial_dimension();

//         has_contact_ = false;
//         if(spatial_dim == 2) {
//             has_contact_ = ConvertContactAlgorithm<2>::apply(params, black_list_, mesh, dof_map, contact_data);
//         } else if(spatial_dim == 3) {
//             has_contact_ = ConvertContactAlgorithm<3>::apply(params, black_list_, mesh, dof_map, contact_data);
//         }

//         if(has_contact_) {
//             contact_tensors_ = contact_data.dof_wise;
//         } else {
//             //init default
//             init_no_contact(
//                 mesh,
//                 dof_map);
//         }

//         has_glue_ = has_contact_ && !params.is_glue->empty();

//         if(has_glue_) {
//             double n_glued = sum(contact_tensors_->is_glue);

//             if(n_glued == 0.) {
//                 has_glue_ = false;
//             }
//         }

//         overall_time.stop();
//         std::cout << "ConvertContactAssembler::assemble: " << overall_time << std::endl;
//         return has_contact_;
//     }

//     void ConvertContactAssembler::couple(const UVector &in, UVector &out) const
//     {
//         assert(contact_tensors_);
//         out = transpose(contact_tensors_->complete_transformation) * in;
//     }

//     void ConvertContactAssembler::uncouple(const UVector &in, UVector &out) const
//     {
//         assert(contact_tensors_);
//         out = contact_tensors_->complete_transformation * in;
//     }

//     void ConvertContactAssembler::couple(const USparseMatrix &in, USparseMatrix &out) const
//     {
//         assert(contact_tensors_);

//         const auto &T = contact_tensors_->complete_transformation;

//         out = transpose(T) * in * T;
//     }

//     const UVector &ConvertContactAssembler::gap() const
//     {
//         assert(contact_tensors_);
//         return contact_tensors_->gap;
//     }

//     UVector &ConvertContactAssembler::gap()
//     {
//         assert(contact_tensors_);
//         return contact_tensors_->gap;
//     }

//     bool ConvertContactAssembler::init_no_contact(
//         const libMesh::MeshBase &mesh,
//         const libMesh::DofMap &dof_map)
//     {

//         if(!contact_tensors_) {
//             contact_tensors_ = std::make_shared<ContactTensors>();
//         }

//         auto n_local_dofs = dof_map.n_local_dofs();

//         contact_tensors_->gap = local_values(n_local_dofs, 100000000);
//         contact_tensors_->weighted_gap = local_values(n_local_dofs, 100000000);
//         contact_tensors_->normal = local_zeros(n_local_dofs);

//         contact_tensors_->inv_mass_vector = local_values(n_local_dofs, 1.);
//         contact_tensors_->is_contact = local_zeros(n_local_dofs);

//         contact_tensors_->T = local_identity(n_local_dofs, n_local_dofs);
//         contact_tensors_->orthogonal_trafo  = local_identity(n_local_dofs, n_local_dofs);
//         contact_tensors_->complete_transformation = local_identity(n_local_dofs, n_local_dofs);
//         has_contact_ = false;
//         has_glue_ = false;
//         return true;
//     }

//     void ConvertContactAssembler::remove_mass(const UVector &in, UVector &out) const
//     {
//         if(!empty(contact_tensors_->inv_mass_vector )) {
//             //P1 stuff
//             if(!empty(contact_tensors_->Q_inv)) {
//                 out = e_mul(
//                     contact_tensors_->inv_mass_vector,
//                     contact_tensors_->Q_inv * in
//                 );

//                 return;
//             }

//             out = e_mul(contact_tensors_->inv_mass_vector, in);

//         } else {
//             out = contact_tensors_->D_inv * in;
//         }
//     }

//     void ConvertContactAssembler::read(Input &in)
//     {
//         in.get("black-list", [this](Input &in) {
//             black_list_ = std::make_shared<ElementBlackList>(true);
//             black_list_->read(in);
//         });
//     }
// }


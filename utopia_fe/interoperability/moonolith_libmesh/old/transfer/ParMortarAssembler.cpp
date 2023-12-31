// #include "ParMortarAssembler.hpp"
// #include "MortarAssemble.hpp"

// #include "libmesh/mesh_inserter_iterator.h"
// #include "libmesh/elem.h"
// #include "libmesh/transient_system.h"
// #include "libmesh/fe.h"

// #include "par_moonolith.hpp"
// #include "moonolith_tree.hpp"
// #include "moonolith_n_tree_mutator_factory.hpp"
// #include "moonolith_n_tree_with_span_mutator_factory.hpp"
// #include "moonolith_n_tree_with_tags_mutator_factory.hpp"
// #include "moonolith_profiler.hpp"
// #include "moonolith_redistribute.hpp"
// #include "moonolith_sparse_matrix.hpp"

// #include <cmath>
// #include <queue>
// #include <algorithm>
// #include <numeric>

// // #include "LibMeshCutlibppAdapters.hpp"
// #include "utopia_ElementDofMap.hpp"
// #include "utopia_VElementAdapter.hpp"
// #include "utopia_VTree.hpp"
// #include "utopia_STree.hpp"

// namespace utopia {

//     using namespace libMesh;

//     ParMortarAssembler::ParMortarAssembler(libMesh::Parallel::Communicator &libmesh_comm,
//                                            const std::shared_ptr<LibMeshFESpaceBase> &master_slave)
//     :libmesh_comm_(libmesh_comm), master_slave_(master_slave)
//     { }

//     class UtopiaMesh {
//     public:

//         explicit UtopiaMesh(const moonolith::Communicator &comm) : comm(comm)
//         {
//             must_destroy_attached[0] = false;
//         }

//         UtopiaMesh(const std::shared_ptr<LibMeshFESpaceBase> &master_slave)
//         {

//             utopiamesh_.reserve(1);

//             utopiamesh_.push_back(master_slave);

//             must_destroy_attached[0] = false;

//             const int n_elements = master_slave->mesh().n_elem();

//             copy_global_dofs(*master_slave, dof_maps_[0], var_type_[0], n_elements, subdomain_id_[0],
//             side_set_id_[0], side_set_id_tag_[0], 101, 102);

//             copy_var_number(*master_slave, var_number_[0]);

//             copy_var_order(*master_slave, var_order_[0]);

//         }

//         inline std::vector< std::shared_ptr<LibMeshFESpaceBase> > &utopiamesh()
//         {
//             return utopiamesh_;

//         }

//         inline const std::vector< std::shared_ptr<LibMeshFESpaceBase> > &utopiamesh() const
//         {
//             return utopiamesh_;

//         }

//         inline long n_elements() const
//         {

//             long ret=0;
//             ret += utopiamesh_[0]->mesh().n_elem();
//             return ret;

//         }

//         inline std::vector<ElementDofMap> &dof_map()
//         {

//             return dof_maps_[0];
//         }

//         inline const std::vector<ElementDofMap> &dof_map() const
//         {

//             return dof_maps_[0];
//         }

//         inline void set_must_destroy_attached(const int index, const bool value)
//         {

//             must_destroy_attached[index] = value;
//         }

//         inline  std::vector<ElementDofMap> &variable_number()
//         {

//             return var_number_[0];
//         }

//         inline const std::vector<ElementDofMap> &variable_number() const
//         {

//             return var_number_[0];
//         }

//         inline std::vector<ElementDofMap> &variable_order()
//         {

//             return var_order_[0];
//         }

//         inline const std::vector<ElementDofMap> &variable_order() const
//         {

//             return var_order_[0];
//         }

//         inline std::vector<ElementDofMap> &variable_type()
//         {

//             return var_type_[0];
//         }

//         inline const std::vector<ElementDofMap> &variable_type() const
//         {
//             return var_type_[0];
//         }

//         inline std::vector<ElementDofMap> &subdomain_id()
//         {

//             return subdomain_id_[0];
//         }

//         inline const std::vector<ElementDofMap> &subdomain_id() const
//         {
//             return subdomain_id_[0];
//         }

//         inline std::vector<ElementDofMap> &side_set_id()
//         {

//             return side_set_id_[0];
//         }

//         inline const std::vector<ElementDofMap> & side_set_id() const
//         {
//             return side_set_id_[0];
//         }

//         inline std::vector<ElementDofMap> &side_set_id_tag()
//         {

//             return side_set_id_tag_[0];
//         }

//         inline const std::vector<ElementDofMap> & side_set_id_tag() const
//         {
//             return side_set_id_tag_[0];
//         }

//     private:
//         moonolith::Communicator comm;
//         std::vector<std::shared_ptr< LibMeshFESpaceBase>> utopiamesh_;
//         std::vector<ElementDofMap> dof_maps_[1];
//         std::vector<ElementDofMap> var_number_[1];
//         std::vector<ElementDofMap> var_order_[1];
//         std::vector<ElementDofMap> var_type_[1];
//         std::vector<ElementDofMap> subdomain_id_[1];
//         std::vector<ElementDofMap> side_set_id_[1];
//         std::vector<ElementDofMap> side_set_id_tag_[1];
//         bool must_destroy_attached[1];

//         inline static void copy_global_dofs(LibMeshFESpaceBase &space, std::vector<ElementDofMap> &dof_map,
//         std::vector<ElementDofMap> &variable_type, const int n_elements,std::vector<ElementDofMap> &subdomain_id,
//         std::vector<ElementDofMap> &side_set_id, std::vector<ElementDofMap> &side_set_id_tag, int tag_1, int tag_2)
//         {

//             auto &mesh = space.mesh();
//             auto &original_dof_map = space.dof_map();
//             std::vector<dof_id_type> temp;

//             MeshBase::const_element_iterator e_it_s = mesh.active_local_elements_begin();
//             const MeshBase::const_element_iterator e_end_s = mesh.active_local_elements_end();

//             dof_map.resize(n_elements);
//             subdomain_id.resize(n_elements);
//             side_set_id.resize(n_elements);
//             side_set_id_tag.resize(n_elements);

//             variable_type.resize(1);

//             bool first=true;

//             std::vector<const Node *> elem_nodes;

//             int jj_side_id_one = 0;
//             int jj_side_id_one_tag = 0;
//              int jj_side_id_one_check = 0;

//             for (; e_it_s != e_end_s; ++e_it_s)
//             {
//                 Elem * elem = *e_it_s;

//                 bool  check_side_id_one=true;
//                 bool  check_side_id_one_tag=true;
//                 bool  check_side_id_one_check=true;

//                 for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){
//                     if (check_side_id_one==true){
//                         side_set_id[elem->id()].global.insert(side_set_id[elem->id()].global.end(),-1);
//                         check_side_id_one=false;
//                         jj_side_id_one++;
//                     }
//                 }

//                 if (elem->on_boundary()){
//                     for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){
//                          if ((mesh.get_boundary_info().has_boundary_id(elem,side_elem,tag_1) ||
//                          mesh.get_boundary_info().has_boundary_id(elem,side_elem,tag_2)) &&
//                          check_side_id_one_tag==true){
//                             side_set_id[elem->id()].global.insert(side_set_id[elem->id()].global.end()-1,mesh.get_boundary_info().boundary_id(elem,side_elem));
//                              check_side_id_one_tag=false;
//                             jj_side_id_one_tag++;
//                         }
//                     }
//                 }

//                 for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){
//                     if (check_side_id_one_check==true){
// //                        std::cout<<"side_set_id[ "<< elem->id() <<" ] =
// "<<side_set_id[elem->id()].global.at(0)<<std::endl;
//                         check_side_id_one_check=false;
//                         jj_side_id_one_check++;
//                     }
//                 }

//                 subdomain_id[elem->id()].global.insert(subdomain_id[elem->id()].global.end(),elem->subdomain_id());
//                 original_dof_map.dof_indices(elem, temp, 0);
// //                std::cout<<" TEMP SIZE "<< temp.size() <<std::endl;
//                 dof_map[elem->id()].global.insert(dof_map[elem->id()].global.end(), temp.begin(), temp.end());

//                 if (first)
//                 {
//                     variable_type[0].global.push_back(elem->type());
//                     first=false;

//                 }

//             }

// //            std::cout<<" jj_side_id_one_check = "<< jj_side_id_one_check <<std::endl;
// //
// //            std::cout<<" jj_side_id_one_tag = "<< jj_side_id_one_tag <<std::endl;

//         }

//         inline static void copy_var_number(LibMeshFESpaceBase &space, std::vector<ElementDofMap> &variable_number)
//         {
//             variable_number.resize(1);
//             variable_number[0].global.push_back(space.var_num());
//         }

//         inline static void copy_var_order(LibMeshFESpaceBase &space, std::vector<ElementDofMap> &variable_order)
//         {
//             variable_order.resize(1);
//             variable_order[0].global.push_back(space.order());
//         }

//     };

//     template<class Iterator>
//     static void write_space(
//         const Iterator &begin,
//         const Iterator &end,
//         LibMeshFESpaceBase &space,
//         const std::vector<ElementDofMap> &dof_map,
//         const std::vector<ElementDofMap> &variable_number,
//         const std::vector<ElementDofMap> &variable_order,
//         const std::vector<ElementDofMap> &subdomain_id,
//         const std::vector<ElementDofMap> &side_set_id,
//         moonolith::OutputStream &os,
//         const int tag_1,
//         const int tag_2)
//     {
//         const int dim 		  = space.mesh().mesh_dimension();
//         const long n_elements = std::distance(begin, end);

//         std::set<long> nodeIds;
//         std::map<long, long> mapping;
//         std::vector<dof_id_type> dof_array;

//         for(Iterator it = begin; it != end; ++it) {

//             const Elem *elem = space.mesh().elem(*it);

//             for(dof_id_type j = 0; j != elem->n_nodes(); ++j) {

//                 nodeIds.insert(elem->node(j));

//             }
//         }

//         long n_nodes = nodeIds.size();

//         // Estimate for allocation
//         os.request_space( (n_elements * 8 + n_nodes * dim) * (sizeof(double) + sizeof(long)) );

//         //WRITE 1
//         os << dim;

//         //std::cout<<"dim write= "<<dim<<std::endl;

//         int index = 0;
//         for (auto nodeId : nodeIds) {
//             mapping[nodeId] = index++;
//         }

//         //WRITE 2
//         os << n_nodes;

//         //WRITE 6
//         os << n_elements;

// //        std::cout<<"elem write= "<<n_elements<<std::endl;

//         for(auto node_id : nodeIds){

//             const Point &p = space.mesh().node(node_id);

//             for(int i = 0; i < dim; ++i) {

//                 //WRITE 3
//                 os << p(i);

//             }

//             //std::cout<<"write_point"<<p<<std::endl;

//         }

//         std::vector<dof_id_type> indices_vector;

//         for(Iterator it = begin; it != end; ++it) {

//             const int k = *it;

//             const Elem *elem = space.mesh().elem(*it);

//             const int e_n_nodes = elem->n_nodes();

//             const int type = elem->type();

//             //WRITE 7
//             os << type << e_n_nodes;

//             for (int i = 0; i != e_n_nodes; ++i) {

//                 auto it = mapping.find(elem->node(i));

//                 assert(it != mapping.end());

//                 int index = it->second;

//                 //WRITE 8
//                 os << index;

//             }

//             //WRITE 9
//             assert(!dof_map.at(elem->id()).empty());

//             os << dof_map.at(elem->id());

//             bool  size=true;

//             int volume_tag;

//             volume_tag=subdomain_id[elem->id()].global.at(0);

//             os << volume_tag;

//             int side_set_tag;

//             bool check_side_id_one=true;

// //            for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){
// //
// //                if(check_side_id_one==true){

//             side_set_tag=side_set_id[elem->id()].global.at(0);

// //                    std::cout<<" write surface role outside= "<< side_set_tag <<std::endl;

//             os << side_set_tag;
// //
// //                    check_side_id_one=false;
// //                }
// //            }
//         }
// //
// //

//         //WRITE 11
//         os << variable_number.at(0);

//         //WRITE 12
//         os << variable_order.at(0);

//     }

//     template<class Iterator>
//     static void write_element_selection(const Iterator &begin, const Iterator &end, const UtopiaMesh &utopiamesh,
//     moonolith::OutputStream &os)
//     {

//         auto m = utopiamesh.utopiamesh()[0];

//         write_space(begin, end, *m, utopiamesh.dof_map(), utopiamesh.variable_number(), utopiamesh.variable_number(),
//         utopiamesh.subdomain_id(), utopiamesh.side_set_id(), os, 101, 102);

//     }

//     static void read_space(moonolith::InputStream &is,
//                            std::shared_ptr<LibMeshFESpaceBase> &space,
//                            std::vector<ElementDofMap> &dof_map,
//                            std::vector<ElementDofMap> &variable_number,
//                            std::vector<ElementDofMap> &variable_order,
//                            std::vector<ElementDofMap> &subdomain_id,
//                            std::vector<ElementDofMap> &side_set_id,
//                            const libMesh::Parallel::Communicator &comm,
//                            int tag_1,
//                            int tag_2)
//     {
//         using namespace std;

//         //READ 1
//         int dim;
//         is >> dim;
//         //std::cout<<"I am reading the mesh "<<std::endl;

//         //READ 2
//         long n_nodes;
//         is >> n_nodes;

//         //READ 6
//         long n_elements;
//         is >> n_elements;

// //        std::cout<<"elem read= "<<n_elements<<std::endl;

//         auto mesh_ptr = std::make_shared<SerialMesh>(comm, dim);

//         mesh_ptr->reserve_nodes(n_nodes);

//         for (long iii = 0; iii != n_nodes; ++iii) {

//             Point p;

//             for(int j = 0; j < dim; ++j) {
//                 //READ 3
//                 is >> p(j);
//             }

//             mesh_ptr->add_point(p);
//             //std::cout<<"read_point"<<p<<std::endl;

//         }

//         dof_map.resize(n_elements);

//         subdomain_id.resize(n_elements);

//         side_set_id.resize(n_elements);

//         for(long i = 0; i !=n_elements; ++i) {

//             //READ 7

//             int type, e_n_nodes;

//             is >> type >> e_n_nodes;

//             //std::cout<<"e_n_nodes_read = "<<e_n_nodes<<std::endl;

//             auto elem =  Elem::build(ElemType(type)).release();

//             //std::cout<<"n_side_read ="<< elem->n_sides()<<std::endl;

//             int index;

//             for (int ii = 0; ii != e_n_nodes; ++ii) {

//                 //READ 8
//                 is >> index;
//                 //std::cout<<"index = "<<index<<std::endl;
//                 elem->set_node(ii) = & mesh_ptr->node(index);

//             }

//             //READ 9
//             is >> dof_map.at(i);
//             //std::cout<< "dof_map_read = "<<dof_map[i].global.at(0)<<std::endl;

//             int volume_tag, side_set_tag;

//             bool on_boundary=false;
//             //std::cout<<"read n_elements = "<<n_elements<<std::endl;

//             is >> volume_tag;

//             //std::cout<<" read volume role = "<< volume_tag <<std::endl;

//             subdomain_id[i].global.insert(subdomain_id[i].global.end(),volume_tag);

// //            bool check_side_id_two=true;
// //
// //            for (unsigned int side=0; side<elem->n_sides(); side++){
// //
// //                if (check_side_id_two==true){

//             is >> side_set_tag;
// //
// //                    if (side_set_tag<0) side_set_tag=1000;
// //
// //                    mesh_ptr->get_boundary_info().add_side(elem, side, side_set_tag);

// //                    std::cout<<" read surface role outside= "<< side_set_tag <<std::endl;
//             side_set_id[i].global.insert(side_set_id[i].global.end(),side_set_tag);

// //                    check_side_id_two=false;
// //                }
// //            }
// //

//             mesh_ptr->add_elem(elem);

//             libmesh_assert(elem);

//         }

//         //READ 11
//         variable_number.resize(1);
//         is >> variable_number.at(0);

//         //READ 12
//         variable_order.resize(1);
//         is >> variable_order.at(0);

//         //EquationSystems es(*mesh_ptr);
//         // const System & new_sys = es.get_system(0);

//         //!!!! dummy parameters
//         space = make_shared<LibMeshFESpaceReal>(mesh_ptr, nullptr, -1);

//     }

//     static void read_spaces(moonolith::InputStream &is, UtopiaMesh &utopiamesh, const libMesh::Parallel::Communicator
//     &comm_mesh)
//     {

//         bool has_master, has_slave;
//         // is >> has_master >> has_slave;

//         utopiamesh.utopiamesh().resize(1);

//         read_space(is, utopiamesh.utopiamesh()[0], utopiamesh.dof_map(), utopiamesh.variable_number(),
//         utopiamesh.variable_order(), utopiamesh.subdomain_id(), utopiamesh.side_set_id(), comm_mesh, 101, 102);

//         utopiamesh.set_must_destroy_attached(0,true);

//     }

//     template<int Dimensions, class Fun>
//     static bool Assemble(moonolith::Communicator &comm,
//                          std::shared_ptr<LibMeshFESpaceBase> &master_slave,
//                          Fun process_fun,
//                          const moonolith::SearchSettings &settings)
//     {
//         using namespace moonolith;

//         typedef VTree<Dimensions> NTreeT;
//         typedef typename NTreeT::DataContainer DataContainer;
//         typedef typename NTreeT::DataType Adapter;

//         const long maxNElements = settings.max_elements;
//         const long maxDepth = settings.max_depth;

//         auto &mesh = master_slave->mesh();

//         const int n_elements = mesh.n_elem();

//         const Parallel::Communicator &libmesh_comm_mesh = master_slave->mesh().comm();

//         MeshBase::const_element_iterator e_it = mesh.active_local_elements_begin();
//         const MeshBase::const_element_iterator e_end = mesh.active_local_elements_end();
//         std::vector<int> block_id;
//         std::vector<int> block_id_def;

//         int i=0;
//         for (; e_it != e_end; ++e_it)
//         {
//             Elem * elem = *e_it;
//             if (i==0){
//                 block_id_def.push_back(elem->subdomain_id());}

//             block_id.push_back(elem->subdomain_id());
//             if (i>0 && block_id.at(i)!=block_id.at(i-1)){
//                 block_id_def.push_back(block_id.at(i));
//             }

//             i++;

//         }

//         auto predicate = std::make_shared<MasterAndSlave>();
//         predicate->add(block_id_def.at(0), block_id_def.at(1));

//         MOONOLITH_EVENT_BEGIN("create_adapters");
//         ////////////////////////////////////////////////////////////////////////////////////////////////////
//         auto tree = NTreeT::New(predicate, maxNElements, maxDepth);
//         tree->reserve(n_elements);

//         auto local_spaces = std::make_shared<UtopiaMesh>(master_slave);

//         for (auto it = master_slave->mesh().active_local_elements_begin();
//         it!=master_slave->mesh().active_local_elements_end(); ++it) {
//             auto elem=*it;
//             int tag = elem->subdomain_id();
//             Adapter a(mesh, elem->id(), elem->id() , tag);
//             assert(!local_spaces->dof_map()[elem->id()].empty());
//             a.set_dof_map(&local_spaces->dof_map()[elem->id()].global);
//             tree->insert(a);
//         }

//         tree->root()->bound().static_bound().enlarge(1e-8);

//         ////////////////////////////////////////////////////////////////////////////////////////////////////
//         MOONOLITH_EVENT_END("create_adapters");

//         //Just to have an indexed-storage
//         std::map<long, std::shared_ptr<UtopiaMesh> > utopiamesh;
//         std::map<long, std::vector<std::shared_ptr<UtopiaMesh> > > migrated_meshes;

//         auto read = [&utopiamesh, &migrated_meshes, block_id, comm, &libmesh_comm_mesh]
//         (
//          const long ownerrank,
//          const long senderrank,
//          bool is_forwarding, DataContainer &data,
//          InputStream &in
//          ) {

//             CHECK_STREAM_READ_BEGIN("vol_proj", in);

//             auto proc_space = std::make_shared<UtopiaMesh>(comm);

//             read_spaces(in, *proc_space, libmesh_comm_mesh);

//             if (!is_forwarding) {
//                 assert(!utopiamesh[ownerrank]);
//                 utopiamesh[ownerrank] = proc_space;
//             } else {
//                 migrated_meshes[ownerrank].push_back(proc_space);
//             }

//             data.reserve(data.size() + proc_space->n_elements());
//             std::vector<long> tag_vec = proc_space->subdomain_id()[0].global;

//             for(auto s : proc_space->utopiamesh()){
//                 for (int i = 0; i<s->mesh().n_elem(); ++i) {
//                     int tag =proc_space->subdomain_id()[i].global.at(0);
//                     data.push_back(Adapter(s->mesh(), i, i, tag));
//                     data.back().set_dof_map(&proc_space->dof_map()[i].global);
//                 }
//             }

//             CHECK_STREAM_READ_END("vol_proj", in);
//         };

//         auto write = [&local_spaces, &utopiamesh, &comm]
//         (
//          const long ownerrank, const long recvrank,
//          const std::vector<long>::const_iterator &begin,
//          const std::vector<long>::const_iterator &end,
//          const DataContainer &data,
//          OutputStream &out) {

//             CHECK_STREAM_WRITE_BEGIN("vol_proj", out);

//             if (ownerrank == comm.rank()) {
//                 write_element_selection(begin, end, *local_spaces, out);

//             } else {
//                 auto it = utopiamesh.find(ownerrank);
//                 assert(it != utopiamesh.end());
//                 std::shared_ptr<UtopiaMesh> spaceptr = it->second;
//                 assert(std::distance(begin, end) > 0);
//                 write_element_selection(begin, end, *spaceptr, out);

//             }

//             CHECK_STREAM_WRITE_END("vol_proj", out);

//         };

//         long n_false_positives = 0, n_intersections = 0;

//         auto fun = [&n_false_positives, &n_intersections, &process_fun](
//                                                                         Adapter &master, Adapter &slave) -> bool {

//             bool ok = process_fun(master, slave);

//             if(ok) {
//                 n_intersections++;
//                 return true;
//             } else {
//                 n_false_positives++;
//                 return false;
//             }
//             return true;

//         };

//         moonolith::search_and_compute(comm, tree, predicate, read, write, fun, settings);

//         long n_total_candidates = n_intersections + n_false_positives;

//         long n_collection[3] = {n_intersections, n_total_candidates, n_false_positives};
//         comm.all_reduce(n_collection, 3, moonolith::MPISum());

//         if (comm.is_root()) {
//             std::cout << "n_intersections: " << n_collection[0]
//             << ", n_total_candidates: " 	 << n_collection[1]
//             << ", n_false_positives: " 	     << n_collection[2] << std::endl;
//         }

//         return true;
//     }

//     template<int Dimensions>
//     bool Assemble(
//                   moonolith::Communicator &comm,
//                   std::shared_ptr<LibMeshFESpaceBase> &master_slave,
//                   USparseMatrix &B,
//                   const moonolith::SearchSettings &settings)
//     {
//         std::shared_ptr<UtopiaMesh> local_fun_spaces = std::make_shared<UtopiaMesh>(master_slave);

//         libMesh::DenseMatrix<libMesh::Real> src_pts;
//         libMesh::DenseMatrix<libMesh::Real> dest_pts;
//         libMesh::DenseMatrix<libMesh::Real> intersection2;
//         Polyhedron src_poly, dest_poly;
//         Polyhedron  intersection3,temp_poly;
//         typedef Intersector::Scalar Scalar;
//         std::shared_ptr<LibMeshFESpaceBase> master_slave_space = master_slave;

//         std::vector<libMesh::dof_id_type> master_dofs, slave_dofs;
//         libMesh::DenseMatrix<libMesh::Real> elemmat;
//         libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;

//         std::shared_ptr<Transform> src_trans;
//         std::shared_ptr<Transform> dest_trans;

//         int skip_zeros = 1;

//         libMesh::Real total_intersection_volume = 0.0;
//         libMesh::Real local_element_matrices_sum = 0.0;

//         //
//         //		moonolith::SparseMatrix<double> mat_buffer(slave->dof_map().n_dofs(),
//         master->dof_map().n_dofs());

//         bool intersected = false;

//         auto fun = [&](const VElementAdapter<Dimensions> &master,
//                        const VElementAdapter<Dimensions> &slave) -> bool {

//             long n_intersections = 0;

//             bool pair_intersected = false;

//             // const auto &src  = master.space();
//             // const auto &dest = slave.space();

//             const auto &src_mesh  = master.space();
//             const auto &dest_mesh = slave.space();

//             // const auto &src_mesh  = src.mesh();
//             // const auto &dest_mesh = dest.mesh();

//             //dest_mesh.print_info();

//             const int src_index  = master.element();
//             const int dest_index = slave.element();

//             auto &src_el  = *src_mesh.elem(src_index);
//             auto &dest_el = *dest_mesh.elem(dest_index);

//             const int dim_src = src_mesh.mesh_dimension();
//             const int dim_sla = dest_mesh.mesh_dimension();

//             if(dim_src == 2)  {
//                 make_polygon(src_el,   src_pts);
//                 make_polygon(dest_el, dest_pts);

//                 if(intersect_2D(src_pts, dest_pts, intersection2)) {

//                     pair_intersected = true;
//                 }
//             }
//             else if(dim_src == 3) {
//                 make_polyhedron(src_el,  src_poly);
//                 make_polyhedron(dest_el, dest_poly);

//                 if(intersect_3D(src_poly, dest_poly, intersection3)) {

//                     pair_intersected = true;

//                 }

//             } else {
//                 assert(false);
//                 return false;
//             }

//             if(pair_intersected) {

//                 return true;

//             } else {

//                 return false;
//             }
//         };

//         if(!Assemble<Dimensions>(comm, master_slave, fun, settings)) {
//             return false;
//         }

//         return true;
//     }

//     bool ParMortarAssembler::Assemble(USparseMatrix &B)
//     {

//         moonolith::SearchSettings settings;

//         moonolith::Communicator comm = libmesh_comm_.get();

//         if(master_slave_->mesh().mesh_dimension() == 2) {
//             //std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
//             return utopia::Assemble<2>(comm, master_slave_, B, settings);
//         }

//         if(master_slave_->mesh().mesh_dimension() == 3) {
//             return utopia::Assemble<3>(comm, master_slave_, B, settings);
//         }

//         assert(false && "Dimension not supported!");
//         return false;
//     }

//     template<int Dimensions, class Fun>
//     static bool SurfaceAssemble(moonolith::Communicator &comm,
//                                 std::shared_ptr<LibMeshFESpaceBase> &master_slave,
//                                 Fun process_fun,
//                                 const moonolith::SearchSettings &settings,
//                                 const libMesh::Real search_radius,
//                                 const int tag_1,
//                                 const int tag_2)
//     {
//         using namespace moonolith;

//         typedef STree<Dimensions> NTreeT;
//         typedef typename NTreeT::DataContainer DataContainer;
//         typedef typename NTreeT::DataType SurfaceAdapter;

//         const long maxNElements = settings.max_elements;
//         const long maxDepth = settings.max_depth;

//         auto &mesh = master_slave->mesh();

//         const int n_elements = mesh.n_elem();
//         const Parallel::Communicator &libmesh_comm_mesh = master_slave->mesh().comm();

//         MeshBase::const_element_iterator e_it = mesh.active_elements_begin();
//         const MeshBase::const_element_iterator e_end = mesh.active_elements_end();
//         std::vector<int> block_id;
//         std::vector<int> block_id_def;

//         int i=0;
//         for (; e_it != e_end; ++e_it)
//         {
//             Elem * elem = *e_it;
//             if (i==0){
//                 block_id_def.push_back(elem->subdomain_id());}

//             block_id.push_back(elem->subdomain_id());
//             if (i>0 && block_id.at(i)!=block_id.at(i-1)){
//                 block_id_def.push_back(block_id.at(i));
//             }

//             i++;
//         }

//         auto predicate = std::make_shared<MasterAndSlave>();
//         predicate->add(tag_1, tag_2);

//         MOONOLITH_EVENT_BEGIN("create_adapters");
//         ////////////////////////////////////////////////////////////////////////////////////////////////////
//         auto tree = NTreeT::New(predicate, maxNElements, maxDepth);
//         tree->reserve(n_elements);

//         auto local_spaces = std::make_shared<UtopiaMesh>(master_slave);

//         int jj=0;

//         for (auto it = mesh.active_local_elements_begin(); it != mesh.active_local_elements_end(); ++it) {
//             auto elem=*it;

//             if(!elem->on_boundary()) {
//                 continue;
//             }

//             bool check_size=false;

//             for(uint side_elem = 0; side_elem < elem->n_sides(); ++side_elem){
//                 if ((predicate->select(master_slave->mesh().get_boundary_info().boundary_id(elem, side_elem))) &&
//                 check_size==false){
//                     SurfaceAdapter a(mesh, elem->id(), elem->id(),
//                     master_slave->mesh().get_boundary_info().boundary_id(elem, side_elem), search_radius);
//                     assert(!local_spaces->dof_map()[elem->id()].empty());
//                     a.set_dof_map(&local_spaces->dof_map()[elem->id()].global);
//                     tree->insert(a);
//                     check_size=true;
//                     jj++;
//                 }
//             }
//         }

//         tree->root()->bound().static_bound().enlarge(1e-8);

//         ////////////////////////////////////////////////////////////////////////////////////////////////////
//         MOONOLITH_EVENT_END("create_adapters");

//         //Just to have an indexed-storage
//         std::map<long, std::shared_ptr<UtopiaMesh> > utopiamesh;
//         std::map<long, std::vector<std::shared_ptr<UtopiaMesh> > > migrated_meshes;

//         auto read = [&utopiamesh, &migrated_meshes, block_id, comm, &libmesh_comm_mesh, search_radius]
//         (
//          const long ownerrank,
//          const long senderrank,
//          bool is_forwarding, DataContainer &data,
//          InputStream &in
//          ) {

//             CHECK_STREAM_READ_BEGIN("vol_proj", in);

//             std::shared_ptr<UtopiaMesh> proc_space = std::make_shared<UtopiaMesh>(comm);

//             read_spaces(in, *proc_space, libmesh_comm_mesh);

//             if (!is_forwarding) {
//                 assert(!utopiamesh[ownerrank]);
//                 utopiamesh[ownerrank] = proc_space;
//             } else {
//                 migrated_meshes[ownerrank].push_back(proc_space);
//             }

//             data.reserve(data.size() + proc_space->n_elements());

//             for(auto s : proc_space->utopiamesh()){
//                 int i=0;
//                 for (int i = 0; i<s->mesh().n_elem(); ++i) {
//                     auto elem=s->mesh().elem(i);
//                     int tag =proc_space->side_set_id()[i].global.at(0);
//                         data.push_back(SurfaceAdapter(s->mesh(), i, i,tag,search_radius));
//                         assert(!proc_space->dof_map()[i].empty());
//                         assert(!proc_space->side_set_id()[i].empty());
//                         data.back().set_dof_map(&proc_space->dof_map()[i].global);
//                 }
//             }

//             CHECK_STREAM_READ_END("vol_proj", in);
//         };

//         auto write = [&local_spaces, &utopiamesh, &comm]
//         (
//          const long ownerrank, const long recvrank,
//          const std::vector<long>::const_iterator &begin,
//          const std::vector<long>::const_iterator &end,
//          const DataContainer &data,
//          OutputStream &out) {

//             CHECK_STREAM_WRITE_BEGIN("vol_proj", out);

//             if (ownerrank == comm.rank()) {
//                 write_element_selection(begin, end, *local_spaces, out);
//             } else {
//                 auto it = utopiamesh.find(ownerrank);
//                 assert(it != utopiamesh.end());
//                 std::shared_ptr<UtopiaMesh> spaceptr = it->second;
//                 assert(std::distance(begin, end) > 0);
//                 write_element_selection(begin, end, *spaceptr, out);

//             }

//             CHECK_STREAM_WRITE_END("vol_proj", out);
//         };

//         long n_false_positives = 0, n_projections = 0;

//         auto fun = [&n_false_positives, &n_projections, &process_fun](SurfaceAdapter &master, SurfaceAdapter &slave)
//         -> bool {
//             bool ok = process_fun(master, slave);

//             if(ok) {
//                 n_projections++;
//                 return true;
//             } else {
//                 n_false_positives++;
//                 return false;
//             }
//             return true;
//         };

//         moonolith::search_and_compute(comm, tree, predicate, read, write, fun, settings);

//         long n_total_candidates = n_projections + n_false_positives;

//         long n_collection[3] = {n_projections, n_total_candidates, n_false_positives};
//         comm.all_reduce(n_collection, 3, moonolith::MPISum());

//         if (comm.is_root()) {
//             std::cout << "n_intersections: " << n_collection[0]
//             << ", n_total_candidates: " 	 << n_collection[1]
//             << ", n_false_positives: " 	     << n_collection[2] << std::endl;
//         }

//         return true;
//     }
//     template<int Dimensions>
//     bool SurfaceAssemble(
//                          moonolith::Communicator &comm,
//                          std::shared_ptr<LibMeshFESpaceBase> &master_slave,
//                          USparseMatrix &B,
//                          const moonolith::SearchSettings &settings,const libMesh::Real search_radius, const int
//                          tag_1, const int tag_2)
//     {

//         assert(false && "never use this!");

//         std::shared_ptr<UtopiaMesh> local_fun_spaces = std::make_shared<UtopiaMesh>(master_slave);

//         libMesh::DenseMatrix<libMesh::Real> src_pts;
//         libMesh::DenseMatrix<libMesh::Real> dest_pts;
//         libMesh::DenseMatrix<libMesh::Real> intersection2;
//         Polyhedron src_poly, dest_poly;
//         Polyhedron  intersection3,temp_poly;

//         std::shared_ptr<LibMeshFESpaceBase> master_slave_space = master_slave;

//         static const double tol = 1e-8;

//         std::vector<libMesh::dof_id_type> master_dofs, slave_dofs;
//         libMesh::DenseMatrix<libMesh::Real> elemmat;
//         libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;
//         DenseMatrix<Real> side_polygon_1, side_polygon_2;
//         DenseMatrix<Real> isect_polygon_1, isect_polygon_2;

//         std::shared_ptr<Transform> src_trans;
//         std::shared_ptr<Transform> dest_trans;

//         Point n1, n2;

//         int skip_zeros = 1;

//         int slave_dof_n=0;

//         int master_dof_n=0;

//         const int dim = master_slave->mesh().mesh_dimension();

//         libMesh::Real total_intersection_volume = 0.0;
//         libMesh::Real local_element_matrices_sum = 0.0;

//         MeshBase::const_element_iterator e_it = master_slave->mesh().active_local_elements_begin();

//         const MeshBase::const_element_iterator e_end = master_slave->mesh().active_local_elements_end();

//         for (; e_it != e_end; ++e_it)
//         {
//             Elem * elem = *e_it;

//             // std::cout <<"subdomain_id = "<< elem->subdomain_id() << "elem_id = " << elem->id() <<std::endl;

//             bool size_new=true;

//             int jj = 0;

//             if (elem->on_boundary()){

//                 for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){

//                     if (master_slave->mesh().get_boundary_info().has_boundary_id(elem,side_elem, tag_1)){

//                         master_dof_n++;
//                     }

//                     if (master_slave->mesh().get_boundary_info().has_boundary_id(elem,side_elem, tag_2)){

//                         slave_dof_n++;
//                     }
//                 }
//             }
//         }

//         int mat_buffer_row = master_slave->dof_map().n_dofs() * master_slave->mesh().n_elem();
//         int mat_buffer_col = master_slave->dof_map().n_dofs() * master_slave->mesh().n_elem();

//         moonolith::SparseMatrix<double> mat_buffer(comm);
//         mat_buffer.set_size( mat_buffer_row, mat_buffer_col);

//         moonolith::SparseMatrix<double> p_buffer(comm);
//         p_buffer.set_size(master_slave->dof_map().n_dofs(), mat_buffer_col);

//         moonolith::SparseMatrix<double> q_buffer(comm);
//         q_buffer.set_size(mat_buffer_row, master_slave->dof_map().n_dofs());

//         moonolith::SparseMatrix<double> rel_area_buff(comm);
//         rel_area_buff.set_size(master_slave->dof_map().n_dofs(), 1);

//         bool intersected = false;
//         auto fun = [&](const SElementAdapter<Dimensions> &master,
//                        const SElementAdapter<Dimensions> &slave) -> bool {

//             long n_intersections = 0;

//             using namespace moonolith;

//             auto predicate = std::make_shared<MasterAndSlave>();
//             predicate->add(tag_1, tag_2);

//             bool pair_intersected = false;

//             //  const auto &src  = master.space();
//             // const auto &dest = slave.space();
//             // const auto &src_mesh  = src.mesh();
//             // const auto &dest_mesh = dest.mesh();

// 			const auto &src_mesh  = master.space();
//             const auto &dest_mesh = slave.space();

//             libMesh::DenseMatrix<libMesh::Real> elemmat;

//             const int src_index  = master.element();
//             const int dest_index = slave.element();

//             auto &src_el  = *src_mesh.elem(src_index);
//             auto &dest_el = *dest_mesh.elem(dest_index);

//             const int dim_src = src_mesh.mesh_dimension();
//             const int dim_sla = dest_mesh.mesh_dimension();

//             Box box_1(dim_src), box_2(dim_sla);

//             QMortar src_ir_ref(dim_src);
//             QMortar src_ir(dim_src);
//             QMortar dest_ir(dim_sla);
//             QMortar dest_ir_ref(dim_sla);

//             std::vector<long> src_order = local_fun_spaces->variable_order()[0].global;
//             const int approx_order=src_order[0];

//             std::shared_ptr<ContactAssembly> surface_assemble;

//             //FIXME This is a hack

//             //
//             //			std::vector<long> src_order = local_fun_spaces->variable_order(0)[0].global;
//             //
//             //			int src_order_new=src_order[0];
//             //
//             //
//             //
//             //			std::vector<long> dest_order = local_fun_spaces->variable_order(1)[0].global;
//             //
//             //			int dest_order_new=dest_order[0];
//             //
//             //
//             //
//             //			std::vector<long> src_var_num = local_fun_spaces->variable_number(0)[0].global;
//             //
//             //			int src_var_new=src_var_num[0];
//             //
//             //
//             //
//             //			std::vector<long> dest_var_num = local_fun_spaces->variable_number(1)[0].global;
//             //
//             //			int dest_var_new=dest_var_num[0];

//             std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
//             //
//             master_fe = libMesh::FEBase::build(src_mesh.mesh_dimension(), FIRST);
//             slave_fe  = libMesh::FEBase::build(dest_mesh.mesh_dimension(), FIRST);
//             //

//             typedef Intersector::Scalar Scalar;

//             if(dim_src == 2)  {
//                 make_polygon(src_el,   src_pts);
//                 make_polygon(dest_el,  dest_pts);
//                 src_trans  = std::make_shared<Transform2>(src_el);
//                 dest_trans = std::make_shared<Transform2>(dest_el);

//             }

//             else if(dim_src == 3) {
//                 make_polyhedron(src_el,  src_poly);
//                 make_polyhedron(dest_el, dest_poly);
//                 src_trans  = std::make_shared<Transform3>(src_el);
//                 dest_trans = std::make_shared<Transform3>(dest_el);

//             }

//             for(uint side_1 = 0; side_1 < src_el.n_sides(); ++side_1) {

//                 libMesh::UniquePtr<libMesh::Elem> s_1 = src_el.side(side_1);

//                 if(!s_1->on_boundary()) continue;

// //                    if(src_el.neighbor_ptr(side_1) != nullptr) continue;
//                 auto side_ptr_1 = src_el.build_side_ptr(side_1);

//                 compute_side_normal(dim_src, *side_ptr_1, n1);

//                 box_1.reset();
//                 enlarge_box_from_side(dim_src, *side_ptr_1, box_1, search_radius);

//                 if(dim_src == 2) {
//                     make_polygon(*side_ptr_1, side_polygon_1);
//                 } else if(dim_src == 3) {
//                     make_polygon_3(*side_ptr_1, side_polygon_1);
//                 } else {
//                     assert(false);
//                 }

//                 for(uint side_2 = 0; side_2 < dest_el.n_sides(); ++side_2) {

//                     //if(dest_el.neighbor_ptr(side_2) != nullptr) continue;

//                     libMesh::UniquePtr<libMesh::Elem> s_2 = dest_el.side(side_2);

//                     if(!s_2->on_boundary()) continue;

//                     auto side_ptr_2 = dest_el.build_side_ptr(side_2);
//                     compute_side_normal(dim_sla, *side_ptr_2, n2);

//                     const Real cos_angle = n1.contract(n2);

//                     //if the angle is more than 60 degrees ( cos(60/180*pi) == 0.5 ) or has same orientation skip
//                     if(cos_angle >= -0.5) {
//                         continue;
//                     }

//                     box_2.reset();
//                     enlarge_box_from_side(dim_sla, *side_ptr_2, box_2, search_radius);

//                     if(!box_1.intersects(box_2, tol)) {
//                         continue;
//                     }

//                     if(dim_sla==2){
//                         make_polygon(*side_ptr_2,side_polygon_2);

//                         if(!project_2D(side_polygon_1,side_polygon_2,isect_polygon_1,isect_polygon_2)){
//                             continue;
//                         }
//                         const Scalar dx = dest_pts(0, 0) - dest_pts(1, 0);
//                         const Scalar dy = dest_pts(0, 1) - dest_pts(1, 1);

//                         const Scalar isect_dx = isect_polygon_2(0, 0) - isect_polygon_2(1, 0);
//                         const Scalar isect_dy = isect_polygon_2(0, 1) - isect_polygon_2(1, 1);

//                         const Scalar area   = std::sqrt(isect_dx*isect_dx + isect_dy*isect_dy);
//                         const Scalar weight = area/std::sqrt(dx*dx + dy*dy);

//                         const int order = order_for_l2_integral(dim_src, src_el, approx_order, dest_el,
//                         approx_order);

//                         make_composite_quadrature_on_surf_2D(isect_polygon_1, weight, order, src_ir);

//                         make_composite_quadrature_on_surf_2D(isect_polygon_2, weight, order, dest_ir);

//                         pair_intersected = true;

//                         surface_assemble = std::make_shared<ContactAssembly>();
//                         surface_assemble->isect_area	   = area;
//                         surface_assemble->relative_area    = weight;

//                     } else if(dim_src == 3) {
//                         make_polygon_3(*side_ptr_2, side_polygon_2);

//                         if(!project_3D(
//                                        side_polygon_1,
//                                        side_polygon_2,
//                                        isect_polygon_1,
//                                        isect_polygon_2))
//                         {
//                             continue;
//                         }

//                         const Scalar area_slave = Intersector::polygon_area_3(side_polygon_2.m(),
//                         &side_polygon_2.get_values()[0]); const Scalar area   	=
//                         Intersector::polygon_area_3(isect_polygon_2.m(), &isect_polygon_2.get_values()[0]); const
//                         Scalar weight 	= area/area_slave;

//                         const int order = order_for_l2_integral(dim_src, src_el, approx_order, dest_el,
//                         approx_order);

//                         make_composite_quadrature_on_surf_3D(isect_polygon_1, weight, order, src_ir);
//                         make_composite_quadrature_on_surf_3D(isect_polygon_2, weight, order, dest_ir);

//                         pair_intersected = true;

//                         surface_assemble = std::make_shared<ContactAssembly>();
//                         surface_assemble->isect_area	= area;
//                         surface_assemble->relative_area = weight;

//                         //use boundary info instead
//                         // rel_area_buff.add(s_2->id(), 0, weight);

//                     } else {
//                         assert(false);
//                         return false;
//                     }
//                 }
//             }

//             if(pair_intersected) {

//                 transform_to_reference_surf(*src_trans,  src_el.type(),  src_ir, src_ir_ref);
//                 transform_to_reference_surf(*dest_trans, dest_el.type(), dest_ir, dest_ir_ref);

//                 master_fe->attach_quadrature_rule(&src_ir_ref);
//                 master_fe->reinit(&src_el);

//                 slave_fe->attach_quadrature_rule(&dest_ir_ref);
//                 slave_fe->reinit(&dest_el);

//                 surface_assemble->parent_element_master  = src_index;

//                 surface_assemble->id_master 			 = src_el.id();

//                 surface_assemble->parent_element_slave   = dest_index;

//                 surface_assemble->id_slave 			     = dest_el.id();

//                 surface_assemble->coupling.zero();

//                 elemmat.zero();

//                 const auto &master_dofs = master.dof_map();
//                 const auto &slave_dofs  = slave.dof_map();

//                 std::vector<dof_id_type> dof_indices_slave_vec(slave_dofs.size());
//                 std::vector<dof_id_type> dof_indices_master_vec(master_dofs.size());

//                 MeshBase::const_element_iterator e_it = master_slave->mesh().active_elements_begin();
//                 const MeshBase::const_element_iterator e_end = master_slave->mesh().active_elements_end();

//                 std::vector<int> block_id;
//                 std::vector<int> block_id_def;

// //
// //                std::cout<<"************************************************** "<<std::endl;
// //
// //
// //                std::cout<<"************* dest_el.id() = "<< dest_el.id() <<std::endl;
// //
// //
// //                std::cout<<"************* src_el.id() = "<< src_el.id() <<std::endl;
// //
// //
// //                std::cout<<"************************************************** "<<std::endl;
// //

//                 int i=0;
//                 for (; e_it != e_end; ++e_it)
//                 {
//                     Elem * elem = *e_it;
//                     if (i==0){
//                         block_id_def.push_back(elem->subdomain_id());}

//                     block_id.push_back(elem->subdomain_id());
//                     if (i>0 && block_id.at(i)!=block_id.at(i-1)){
//                         block_id_def.push_back(block_id.at(i));
//                     }

//                     i++;

//                 }

//                 for(uint i = 0; i <  slave_dofs.size(); ++i) {

//                     Elem *elem_off = *( master_slave->mesh().active_subdomain_elements_begin(block_id_def.at(1)));

// //                    std::cout<<"************* dest_el.id() = "<< dest_el.id() <<std::endl;
// //
// //                    std::cout<<"************* elem_off->id() = "<< elem_off->id()<<std::endl;

//                     dof_indices_slave_vec[i] =  dest_el.id() * dest_el.n_nodes() + i;

// //                    std::cout<< "************ dof_I = "<<  dest_el.id() * dest_el.n_nodes() + i << std::endl;
// //
// //                    std::cout<< "************ dof_I_comparison = "<<  slave_dofs[i] << std::endl;
//                 }

//                 for(uint i = 0; i <  master_dofs.size(); ++i) {

//                     Elem *elem_off = *( master_slave->mesh().active_subdomain_elements_begin(block_id_def.at(0)));

// //                    std::cout<<"************* src_el.id() = "<< src_el.id() <<std::endl;

// //                    std::cout<<"************* elem_off->id() = "<< elem_off->id()<<std::endl;

//                     dof_indices_master_vec[i] = src_el.id() * src_el.n_nodes() + i;

// //                    std::cout<< "************ dof_J = "<< src_el.id() * src_el.n_nodes() + i << std::endl;
// //
// //                    std::cout<< "************ dof_J_comparison = "<<  master_dofs[i] << std::endl;

//                 }

//                 mortar_assemble(*master_fe, *slave_fe, elemmat);

//                 // std::cout << "-----------------------------------------\n";
//                 // std::cout << src_index << ", " << dest_index << "\n";
//                 // elemmat.print(std::cout);
//                 // for(auto i : slave_dofs) {
//                 // 	std::cout << i << " ";
//                 // }
//                 // std::cout << "\n";

//                 // for(auto i : master_dofs) {
//                 // 	std::cout << i << " ";
//                 // }
//                 // std::cout << "\n";
//                 // std::cout << "-----------------------------------------\n";

//                auto partial_sum = std::accumulate(elemmat.get_values().begin(), elemmat.get_values().end(),
//                libMesh::Real(0.0));

//                // const Scalar local_mat_sum = std::accumulate(surface_assemble->coupling.get_values().begin(),
//                surface_assemble->coupling.get_values().end(), libMesh::Real(0.0));

//                local_element_matrices_sum += partial_sum;

//                assert(slave_dofs.size() == elemmat.m());
//                assert(master_dofs.size() == elemmat.n());

//                 // std::cout << src_index << ", " << dest_index << ": " << partial_sum << std::endl;
//                 // dest_ir.print_info();

//                 //std::cout << "elemmat.n =" << elemmat.n() << "\n";

//                 //std::cout << "elemmat.m =" << elemmat.m() << "\n";

//                 for(int i = 0; i <  slave_dofs.size(); ++i) {
//                     const long dof_I = slave_dofs[i];
//                     const long dof_J = dof_indices_slave_vec[i];
//                     p_buffer.set(dof_I, dof_J, 1.);
//                 }

//                 for(int i = 0; i <  dof_indices_master_vec.size(); ++i) {
//                     const long dof_I = dof_indices_master_vec[i];
//                     const long dof_J = master_dofs[i];
//                     q_buffer.set(dof_I, dof_J, 1.);
//                 }

//                 for(int i = 0; i <  dof_indices_slave_vec.size(); ++i) {

//                     const long dof_I = dof_indices_slave_vec[i];

//                     //std::cout<< "************ dof_I_index = "<< dof_I <<std::endl;

//                     for(int j = 0; j <  dof_indices_master_vec.size(); ++j) {

//                         const long dof_J = dof_indices_master_vec[j];

//                         //std::cout<< "************ dof_J_index = "<< dof_J <<std::endl;

//                         mat_buffer.add(dof_I, dof_J, elemmat(i, j));
//                     }
//                 }

//                 int dim = dim_src;

// //                std::cout<< "surface_assemble->isect_area = " << surface_assemble->isect_area <<std::endl;
// //
// //                std::cout<<" pow(surface_assemble->isect_area, dim/(dim-1.)) * dim = " <<
// pow(surface_assemble->isect_area, dim/(dim-1.)) * dim  <<std::endl;
// //
// //                std::cout<<" partial sum = " << partial_sum <<std::endl;
// //
// //                std::cout<<" local_element_matrices_sum = " << local_element_matrices_sum <<std::endl;
// //
// //

//                 assert(fabs(partial_sum - pow(surface_assemble->isect_area, dim/(dim-1.))) < 1e-8 ||
//                 (!is_quad(dest_el.type()) && !is_hex(dest_el.type())));

//                 return true;

//                 //Assemble p
//                 //Iteri su elementi
//                 //calcoli i=dof_map.index(elemen), j=el_id * n_dofs_x_el + local_dof_offset (local_dof_id e.g.,
//                 triangolo = {0, 1, 2}, quad = {0, 1, 2, 3})

//             } else {

//                 return false;
//             }

//         };

//         if(!SurfaceAssemble<Dimensions>(comm, master_slave, fun, settings, search_radius, tag_1, tag_2)) {
//             return false;
//         }

//         // std::cout << mat_buffer << std::endl;

//         double volumes[1] = { local_element_matrices_sum };

//         comm.all_reduce(volumes, 1, moonolith::MPISum());

//         const processor_id_type master_proc_id  = master_slave->mesh().processor_id();

//         const dof_id_type n_dofs_on_proc_master = master_slave->dof_map().n_dofs_on_processor(master_proc_id) *
//         master_slave->mesh().n_local_elem();

//         const processor_id_type slave_proc_id   = master_slave->mesh().processor_id();

//         const dof_id_type n_dofs_on_proc_slave  = master_slave->dof_map().n_dofs_on_processor(slave_proc_id) *
//         master_slave->mesh().n_local_elem();

//         if(comm.is_root()) {
//             std::cout << "sum(B): " << volumes[0] <<std::endl;
//         }

//         std::vector<moonolith::Integer>  ownershipRangesMaster(comm.size()+1, 0);
//         std::vector<moonolith::Integer>  ownershipRangesSlave(comm.size()+1, 0);

//         ownershipRangesMaster[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_master);
//         ownershipRangesSlave[comm.rank()+1]  += static_cast<unsigned int>(n_dofs_on_proc_slave);

//         comm.all_reduce(&ownershipRangesMaster[0], ownershipRangesMaster.size(), moonolith::MPISum());
//         comm.all_reduce(&ownershipRangesSlave[0],  ownershipRangesSlave.size(),  moonolith::MPISum());

//         std::partial_sum(ownershipRangesMaster.begin(), ownershipRangesMaster.end(),
//                          ownershipRangesMaster.begin());

//         std::partial_sum(ownershipRangesSlave.begin(), ownershipRangesSlave.end(),
//                          ownershipRangesSlave.begin());

//         moonolith::Redistribute< moonolith::SparseMatrix<double> > redist(comm.get_mpi_comm());

//         redist.apply(ownershipRangesSlave, mat_buffer, moonolith::AddAssign<double>());

//         assert(ownershipRangesSlave.empty() == ownershipRangesMaster.empty() || ownershipRangesMaster.empty());

//         moonolith::root_describe("petsc assembly begin", comm, std::cout);

//         SizeType  mMaxRowEntries = mat_buffer.local_max_entries_x_col();

//         comm.all_reduce(&mMaxRowEntries, 1, moonolith::MPIMax());

//         const SizeType local_range_slave_range  = ownershipRangesSlave [comm.rank()+1] - ownershipRangesSlave
//         [comm.rank()]; const SizeType local_range_master_range = ownershipRangesMaster[comm.rank()+1] -
//         ownershipRangesMaster[comm.rank()];

//         USparseMatrix B_tilde = utopia::local_sparse(local_range_slave_range, local_range_master_range,
//         mMaxRowEntries);

//         UVector relAreaVec = zeros(master_slave->dof_map().n_dofs());

//         {
//             utopia::Write<utopia::UVector> write(relAreaVec);
//             for (auto it = rel_area_buff.iter(); it; ++it) {
//                 relAreaVec.add(it.row(), *it);
//             }
//         }

//         // disp(relAreaVec);

//         {
//             utopia::Write<utopia::USparseMatrix> write(B_tilde);
//             for (auto it = mat_buffer.iter(); it; ++it) {
//                 B_tilde.set(it.row(), it.col(), *it);

//             }
//         }

//         /*P_tilde*/

//         const dof_id_type n_dofs_on_proc_master_p = master_slave->dof_map().n_dofs_on_processor(master_proc_id) *
//         master_slave->mesh().n_local_elem();

//         const dof_id_type n_dofs_on_proc_slave_p  = master_slave->dof_map().n_dofs_on_processor(slave_proc_id);

//         if(comm.is_root()) {
//             std::cout << "sum(B): " << volumes[0] <<std::endl;
//         }

//         std::vector<moonolith::Integer>  ownershipRangesMaster_p(comm.size()+1, 0);
//         std::vector<moonolith::Integer>  ownershipRangesSlave_p(comm.size()+1, 0);

//         ownershipRangesMaster_p[comm.rank()+1]+= static_cast<unsigned int>(n_dofs_on_proc_master_p);
//         ownershipRangesSlave_p[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_slave_p);

//         comm.all_reduce(&ownershipRangesMaster_p[0], ownershipRangesMaster_p.size(), moonolith::MPISum());

//         comm.all_reduce(&ownershipRangesSlave_p[0],  ownershipRangesSlave_p.size(),  moonolith::MPISum());

//         std::partial_sum(ownershipRangesMaster_p.begin(), ownershipRangesMaster_p.end(),
//                          ownershipRangesMaster_p.begin());

//         std::partial_sum(ownershipRangesSlave_p.begin(), ownershipRangesSlave_p.end(),
//                          ownershipRangesSlave_p.begin());

//         redist.apply(ownershipRangesSlave_p, p_buffer, moonolith::AddAssign<double>());

//         assert(ownershipRangesSlave_p.empty() == ownershipRangesMaster_p.empty() || ownershipRangesMaster_p.empty());

//         moonolith::root_describe("petsc assembly begin", comm, std::cout);

//         SizeType  mMaxRowEntries_p = p_buffer.local_max_entries_x_col();

//         comm.all_reduce(&mMaxRowEntries_p, 1, moonolith::MPIMax());

//         const SizeType local_range_slave_range_p  = ownershipRangesSlave_p [comm.rank()+1] - ownershipRangesSlave_p
//         [comm.rank()]; const SizeType local_range_master_range_p = ownershipRangesMaster_p[comm.rank()+1] -
//         ownershipRangesMaster_p[comm.rank()];

//         USparseMatrix P_tilde = utopia::local_sparse(local_range_slave_range_p, local_range_master_range_p,
//         mMaxRowEntries_p);

//         {
//             utopia::Write<utopia::USparseMatrix> write(P_tilde);
//             for (auto it = p_buffer.iter(); it; ++it) {
//                 P_tilde.set(it.row(), it.col(), *it);

//             }
//         }

//         /*Q_transpose*/

//         const dof_id_type n_dofs_on_proc_master_q = master_slave->dof_map().n_dofs_on_processor(master_proc_id);

//         const dof_id_type n_dofs_on_proc_slave_q  = master_slave->dof_map().n_dofs_on_processor(slave_proc_id) *
//         master_slave->mesh().n_local_elem();

//         if(comm.is_root()) {
//             std::cout << "sum(B): " << volumes[0] <<std::endl;
//         }

//         std::vector<moonolith::Integer>  ownershipRangesMaster_q(comm.size()+1, 0);
//         std::vector<moonolith::Integer>  ownershipRangesSlave_q(comm.size()+1, 0);

//         ownershipRangesMaster_q[comm.rank()+1]+= static_cast<unsigned int>(n_dofs_on_proc_master_q);

//         ownershipRangesSlave_q[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_slave_q);

//         comm.all_reduce(&ownershipRangesMaster_q[0], ownershipRangesMaster_q.size(), moonolith::MPISum());

//         comm.all_reduce(&ownershipRangesSlave_q[0],  ownershipRangesSlave_q.size(),  moonolith::MPISum());

//         std::partial_sum(ownershipRangesMaster_q.begin(), ownershipRangesMaster_q.end(),
//                          ownershipRangesMaster_q.begin());

//         std::partial_sum(ownershipRangesSlave_q.begin(), ownershipRangesSlave_q.end(),
//                          ownershipRangesSlave_q.begin());

//         redist.apply(ownershipRangesSlave_q, q_buffer, moonolith::AddAssign<double>());

//         assert(ownershipRangesSlave_q.empty() == ownershipRangesMaster_q.empty() || ownershipRangesMaster_q.empty());

//         moonolith::root_describe("petsc assembly begin", comm, std::cout);

//         SizeType  mMaxRowEntries_q = q_buffer.local_max_entries_x_col();

//         comm.all_reduce(&mMaxRowEntries_q, 1, moonolith::MPIMax());

//         const SizeType local_range_slave_range_q  = ownershipRangesSlave_q [comm.rank()+1] - ownershipRangesSlave_q
//         [comm.rank()]; const SizeType local_range_master_range_q = ownershipRangesMaster_q[comm.rank()+1] -
//         ownershipRangesMaster_q[comm.rank()];

//         USparseMatrix Q_transpose = utopia::local_sparse(local_range_slave_range_q, local_range_master_range_q,
//         mMaxRowEntries);

//         {
//             utopia::Write<utopia::USparseMatrix> write(Q_transpose);
//             for (auto it = q_buffer.iter(); it; ++it) {
//                 Q_transpose.set(it.row(), it.col(), *it);

//             }
//         }

//        B = P_tilde * B_tilde * Q_transpose;

//        moonolith::root_describe("petsc assembly end", comm, std::cout);
//       return true;
//     }

//     bool ParMortarAssembler::SurfaceAssemble(USparseMatrix &B, const libMesh::Real search_radius, const int tag_1,
//     const int tag_2)
//     {

//         moonolith::SearchSettings settings;

//         moonolith::Communicator comm = libmesh_comm_.get();

//         if(master_slave_->mesh().mesh_dimension() == 2) {
//             //std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
//             return utopia::SurfaceAssemble<2>(comm, master_slave_, B, settings, search_radius, tag_1, tag_2);
//         }

//         if(master_slave_->mesh().mesh_dimension() == 3) {
//             return utopia::SurfaceAssemble<3>(comm, master_slave_, B, settings, search_radius, tag_1, tag_2);
//         }

//         assert(false && "Dimension not supported!");
//         return false;
//     }

//     //
//     //
//     //
//     //	bool ParMortarAssembler::Transfer(USparseMatrix &B, USparseMatrix &T)
//     //	{
//     //		UVector               diag_elem;
//     //		PetscInt               mG;
//     //		PetscInt               nG;
//     //		PetscInt               mL;
//     //		PetscInt               nL;
//     //
//     //		diag_elem = 1./sum(B,1);
//     //		T = diag(diag_elem)*B;
//     //
//     //		return true;
//     //	}
//     //
// }

// //    template<int Dimensions>
// //    bool SurfaceAssemble(
// //                         moonolith::Communicator &comm,
// //                         std::shared_ptr<LibMeshFESpaceBase> &master_slave,
// //                         USparseMatrix &B,
// //                         const moonolith::SearchSettings &settings,const libMesh::Real search_radius, const int
// tag_1, const int tag_2)
// //    {
// //        std::shared_ptr<UtopiaMesh> local_fun_spaces = std::make_shared<UtopiaMesh>(master_slave);
// //
// //        libMesh::DenseMatrix<libMesh::Real> src_pts;
// //        libMesh::DenseMatrix<libMesh::Real> dest_pts;
// //        libMesh::DenseMatrix<libMesh::Real> intersection2;
// //        Polyhedron src_poly, dest_poly;
// //        Polyhedron  intersection3,temp_poly;
// //
// //        std::shared_ptr<LibMeshFESpaceBase> master_slave_space = master_slave;
// //
// //        static const double tol = 1e-8;
// //
// //
// //        std::vector<libMesh::dof_id_type> master_dofs, slave_dofs;
// //        libMesh::DenseMatrix<libMesh::Real> elemmat;
// //        libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;
// //        DenseMatrix<Real> side_polygon_1, side_polygon_2;
// //        DenseMatrix<Real> isect_polygon_1, isect_polygon_2;
// //
// //        std::shared_ptr<Transform> src_trans;
// //        std::shared_ptr<Transform> dest_trans;
// //
// //        Point n1, n2;
// //
// //
// //
// //        int skip_zeros = 1;
// //
// //
// //        libMesh::Real total_intersection_volume = 0.0;
// //        libMesh::Real local_element_matrices_sum = 0.0;
// //
// //        //
// //        //		moonolith::SparseMatrix<double> mat_buffer(slave->dof_map().n_dofs(),
// master->dof_map().n_dofs());
// //
// //        bool intersected = false;
// //
// //        auto fun = [&](const SElementAdapter<Dimensions> &master,
// //                       const SElementAdapter<Dimensions> &slave) -> bool {
// //
// //
// //            using namespace moonolith;
// //            using namespace moonolith;
// //            using namespace moonolith;
// //
// //            auto predicate = make_shared<MasterAndSlave>();
// //            predicate->add(tag_1,tag_2);
// //
// //            long n_intersections = 0;
// //
// //            //std::cout<<"ciao sn in fun"<<std::endl;
// //
// //            bool pair_intersected = false;
// //
// //            const auto &src  = master.space();
// //            const auto &dest = slave.space();
// //
// //            const auto &src_mesh  = src.mesh();
// //            const auto &dest_mesh = dest.mesh();
// //
// //            //dest_mesh.print_info();
// //
// //            const int src_index  = master.element();
// //            const int dest_index = slave.element();
// //
// //            auto &src_el  = *src_mesh.elem(src_index);
// //            auto &dest_el = *dest_mesh.elem(dest_index);
// //
// //            const int dim_src = src_mesh.mesh_dimension();
// //            const int dim_sla = dest_mesh.mesh_dimension();
// //            Box box_1(dim_src), box_2(dim_sla);
// //
// //            QMortar src_ir_ref(dim_src);
// //            QMortar src_ir(dim_src);
// //            QMortar dest_ir(dim_sla);
// //            QMortar dest_ir_ref(dim_sla);
// //
// //
// //
// //
// //            std::vector<long> src_order = local_fun_spaces->variable_order()[0].global;
// //            const int approx_order=src_order[0];
// //
// //
// //            //FIXME This is a hack
// //
// //            //
// //            //			std::vector<long> src_order = local_fun_spaces->variable_order(0)[0].global;
// //            //
// //            //			int src_order_new=src_order[0];
// //            //
// //            //
// //            //
// //            //			std::vector<long> dest_order = local_fun_spaces->variable_order(1)[0].global;
// //            //
// //            //			int dest_order_new=dest_order[0];
// //            //
// //            //
// //            //
// //            //			std::vector<long> src_var_num = local_fun_spaces->variable_number(0)[0].global;
// //            //
// //            //			int src_var_new=src_var_num[0];
// //            //
// //            //
// //            //
// //            //			std::vector<long> dest_var_num = local_fun_spaces->variable_number(1)[0].global;
// //            //
// //            //			int dest_var_new=dest_var_num[0];
// //
// //
// //
// //
// //            std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
// //            //
// //            master_fe = libMesh::FEBase::build(src_mesh.mesh_dimension(), FIRST);
// //            slave_fe  = libMesh::FEBase::build(dest_mesh.mesh_dimension(), FIRST);
// //            //
// //
// //            typedef Intersector::Scalar Scalar;
// //            //
// //            //
// //            //			const int order = order_for_l2_integral(dim, src_el, src_order_new, dest_el,
// dest_order_new);
// //
// //            //            if(!has_constrained_dofs(src, src_el) &&
// //            //               !has_constrained_dofs(dest, dest_el))
// //            {
// //
// //                if(dim_src == 2)  {
// //                    make_polygon(src_el,   src_pts);
// //                    make_polygon(dest_el,  dest_pts);
// //                    {
// //                        //					total_intersection_volume +=
// fabs(Intersector::polygon_area_2(intersection2.m(), &intersection2.get_values()[0]));
// //                        //
// //                        //					const libMesh::Real
// weight=Intersector::polygon_area_2(dest_pts.m(), &dest_pts.get_values()[0]);
// //                        //
// //                        //					make_composite_quadrature_2D(intersection2, weight, order,
// composite_ir);
// //                        //					pair_intersected = true;
// //                        //
// //                        src_trans  = std::make_shared<Transform2>(src_el);
// //                        dest_trans = std::make_shared<Transform2>(dest_el);
// //                        //pair_intersected = true;
// //                    }
// //                }
// //
// //                else if(dim_src == 3) {
// //                    make_polyhedron(src_el,  src_poly);
// //                    make_polyhedron(dest_el, dest_poly);
// //                    {
// //                        //					total_intersection_volume +=
// Intersector::p_mesh_volume_3(intersection3);
// //                        //
// //                        //					const libMesh::Real weight =
// Intersector::p_mesh_volume_3(dest_poly);
// //                        //
// //                        //
// //                        //					make_composite_quadrature_3D(intersection3, weight, order,
// composite_ir);
// //                        src_trans  = std::make_shared<Transform3>(src_el);
// //                        dest_trans = std::make_shared<Transform3>(dest_el);
// //                        //pair_intersected = true;
// //                    }
// //                }
// //
// //                for(uint side_1 = 0; side_1 < src_el.n_sides(); ++side_1) {
// //
// //                    if(src_el.neighbor_ptr(side_1) != nullptr) continue;
// //                    auto side_ptr_1 = src_el.build_side_ptr(side_1);
// //
// //                    compute_side_normal(dim_src, *side_ptr_1, n1);
// //
// //                    box_1.reset();
// //                    enlarge_box_from_side(dim_src, *side_ptr_1, box_1, search_radius);
// //
// //                    if(dim_src == 2) {
// //                        make_polygon(*side_ptr_1, side_polygon_1);
// //                    } else if(dim_src == 3) {
// //                        make_polygon_3(*side_ptr_1, side_polygon_1);
// //                    } else {
// //                        assert(false);
// //                    }
// //
// //
// //
// //                    for(uint side_2 = 0; side_2 < dest_el.n_sides(); ++side_2) {
// //
// //                        if(dest_el.neighbor_ptr(side_2) != nullptr) continue;
// //
// //                        auto side_ptr_2 = dest_el.build_side_ptr(side_2);
// //                        compute_side_normal(dim_sla, *side_ptr_2, n2);
// //
// //                        const Real cos_angle = n1.contract(n2);
// //
// //                        //if the angle is more than 60 degrees ( cos(60/180*pi) == 0.5 ) or has same orientation
// skip
// //                        if(cos_angle >= -0.5) {
// //                            continue;
// //                        }
// //
// //                        box_2.reset();
// //                        enlarge_box_from_side(dim_sla, *side_ptr_2, box_2, search_radius);
// //
// //                        if(!box_1.intersects(box_2, tol)) {
// //                            continue;
// //                        }
// //
// //                        if(dim_sla==2){
// //                            make_polygon(*side_ptr_2,side_polygon_2);
// //                            if(!project_2D(side_polygon_1,side_polygon_2,isect_polygon_1,isect_polygon_2)){
// //                                continue;
// //                            }
// //                            const Scalar dx = dest_pts(0, 0) - dest_pts(1, 0);
// //                            const Scalar dy = dest_pts(0, 1) - dest_pts(1, 1);
// //
// //                            const Scalar isect_dx = isect_polygon_2(0, 0) - isect_polygon_2(1, 0);
// //                            const Scalar isect_dy = isect_polygon_2(0, 1) - isect_polygon_2(1, 1);
// //
// //                            const Scalar area   = std::sqrt(isect_dx*isect_dx + isect_dy*isect_dy);
// //                            const Scalar weight = area/std::sqrt(dx*dx + dy*dy);
// //
// //                            const int order = order_for_l2_integral(dim_src, src_el, approx_order, dest_el,
// approx_order);
// //
// //                            make_composite_quadrature_on_surf_2D(isect_polygon_1, weight, order, src_ir);
// //                            make_composite_quadrature_on_surf_2D(isect_polygon_2, weight, order, dest_ir);
// //
// //                            pair_intersected = true;
// //
// //
// //                            ++n_intersections;
// //
// //
// //                            //                            current_contact = std::make_shared<ContactAssembly>();
// //                            //                            current_contact->isect_area	   = area;
// //                            //                            current_contact->relative_area = weight;
// //                        } else if(dim_src == 3) {
// //                            make_polygon_3(*side_ptr_2, side_polygon_2);
// //
// //                            if(!project_3D(
// //                                           side_polygon_1,
// //                                           side_polygon_2,
// //                                           isect_polygon_1,
// //                                           isect_polygon_2))
// //                            {
// //                                continue;
// //                            }
// //
// //                            const Scalar area_slave = Intersector::polygon_area_3(side_polygon_2.m(),
// &side_polygon_2.get_values()[0]);
// //                            const Scalar area   	= Intersector::polygon_area_3(isect_polygon_2.m(),
// &isect_polygon_2.get_values()[0]);
// //                            const Scalar weight 	= area/area_slave;
// //
// //                            const int order = order_for_l2_integral(dim_src, src_el, approx_order, dest_el,
// approx_order);
// //
// //                            make_composite_quadrature_on_surf_3D(isect_polygon_1, weight, order, src_ir);
// //                            make_composite_quadrature_on_surf_3D(isect_polygon_2, weight, order, dest_ir);
// //
// //                            pair_intersected = true;
// //
// //
// //                            ++n_intersections;
// //
// //                            //
// //                            //                            current_contact = std::make_shared<ContactAssembly>();
// //                            //                            current_contact->isect_area	   = area;
// //                            //                            current_contact->relative_area = weight;
// //                        } else {
// //                            assert(false);
// //                            return false;
// //                        }
// //                    }
// //                }
// //            }
// //
// //
// //
// //
// //
// //
// //            //                   bool pair_intersected = false;
// //            //
// //            //                   if(dim == 2) {
// //            //                            make_polygon(*side_ptr_2, side_polygon_2);
// //            //
// //            //                            if(!project_2D(side_polygon_1, side_polygon_2, isect_polygon_1,
// isect_polygon_2)) {
// //            //                                continue;
// //            //                            }
// //            //
// //            //                            const Scalar dx = polygon_2(0, 0) - polygon_2(1, 0);
// //            //                            const Scalar dy = polygon_2(0, 1) - polygon_2(1, 1);
// //            //
// //            //                            const Scalar isect_dx = isect_polygon_2(0, 0) - isect_polygon_2(1, 0);
// //            //                            const Scalar isect_dy = isect_polygon_2(0, 1) - isect_polygon_2(1, 1);
// //            //
// //            //                            const Scalar area   = std::sqrt(isect_dx*isect_dx + isect_dy*isect_dy);
// //            //                            const Scalar weight = area/std::sqrt(dx*dx + dy*dy);
// //            //
// //            //                            const int order = order_for_l2_integral(dim, el_1, approx_order, el_2,
// approx_order);
// //            //
// //            //                            make_composite_quadrature_on_surf_2D(isect_polygon_1, weight, order,
// q_1);
// //            //                            make_composite_quadrature_on_surf_2D(isect_polygon_2, weight, order,
// q_2);
// //            //
// //            //                            pair_intersected = true;
// //            //                            ++n_projections;
// //            //
// //            //
// //            //                            current_contact = std::make_shared<ContactAssembly>();
// //            //                            current_contact->isect_area	   = area;
// //            //                            current_contact->relative_area = weight;
// //            //
// //            //
// //            //                        } else if(dim == 3) {
// //            //                            make_polygon_3(*side_ptr_2, side_polygon_2);
// //            //
// //            //                            if(!project_3D(
// //            //                                           side_polygon_1,
// //            //                                           side_polygon_2,
// //            //                                           isect_polygon_1,
// //            //                                           isect_polygon_2))
// //            //                            {
// //            //                                continue;
// //            //                            }
// //            //
// //            //                            const Scalar area_slave = Intersector::polygon_area_3(side_polygon_2.m(),
// &side_polygon_2.get_values()[0]);
// //            //                            const Scalar area   	=
// Intersector::polygon_area_3(isect_polygon_2.m(), &isect_polygon_2.get_values()[0]);
// //            //                            const Scalar weight 	= area/area_slave;
// //            //
// //            //                            const int order = order_for_l2_integral(dim, el_1, approx_order, el_2,
// approx_order);
// //            //
// //            //                            make_composite_quadrature_on_surf_3D(isect_polygon_1, weight, order,
// q_1);
// //            //                            make_composite_quadrature_on_surf_3D(isect_polygon_2, weight, order,
// q_2);
// //            //
// //            //                            pair_intersected = true;
// //            //                            ++n_projections;
// //            //
// //            //
// //            //                            current_contact = std::make_shared<ContactAssembly>();
// //            //                            current_contact->isect_area	   = area;
// //            //                            current_contact->relative_area = weight;
// //
// //            //                        } else {
// //            //                            assert(false);
// //            //                            return false;
// //            //                        }
// //
// //
// //            if(pair_intersected) {
// //
// //
// ////                transform_to_reference_surf(*src_trans,  src_el.type(),  src_ir, src_ir_ref);
// ////                transform_to_reference_surf(*dest_trans, dest_el.type(), dest_ir, dest_ir_ref);
// //                //
// //                //               // src.dof_map().dof_indices(&src_el,  master_dofs);
// //                //               // dest.dof_map().dof_indices(&dest_el, slave_dofs);
// //                //
// //                //
// //                //				const auto &master_dofs = master.dof_map();
// //                //				const auto &slave_dofs  = slave.dof_map();
// //                //               // composite_ir.print_info();
// //                //
// //                //
// ////                master_fe->attach_quadrature_rule(&src_ir_ref);
// ////                master_fe->reinit(&src_el);
// ////
// ////                slave_fe->attach_quadrature_rule(&dest_ir_ref);
// ////                slave_fe->reinit(&dest_el);
// //                //
// //                //				elemmat.zero();
// //                //
// //                //				mortar_assemble(*master_fe, *slave_fe, elemmat);
// //                //
// //                //				// std::cout << "-----------------------------------------\n";
// //                //				// std::cout << src_index << ", " << dest_index << "\n";
// //                //				// elemmat.print(std::cout);
// //                //				// for(auto i : slave_dofs) {
// //                //				// 	std::cout << i << " ";
// //                //				// }
// //                //				// std::cout << "\n";
// //                //
// //                //				// for(auto i : master_dofs) {
// //                //				// 	std::cout << i << " ";
// //                //				// }
// //                //				// std::cout << "\n";
// //                //				// std::cout << "-----------------------------------------\n";
// //                //
// //                //				auto partial_sum = std::accumulate(elemmat.get_values().begin(),
// elemmat.get_values().end(), libMesh::Real(0.0));
// //                //                // std::cout << src_index << ", " << dest_index << ": " << partial_sum <<
// std::endl;
// //                //                // dest_ir.print_info();
// //                //
// //                //				local_element_matrices_sum += partial_sum;
// //                //
// //                //				intersected = true;
// //                //
// //                //				++n_intersections;
// //                //std::cout<<"n_intersections"<<n_intersections<<std::endl;
// //                //
// //                //				assert(slave_dofs.size() == elemmat.m());
// //                //				assert(master_dofs.size() == elemmat.n());
// //                //
// //                //				for(int i = 0; i < slave_dofs.size(); ++i) {
// //                //
// //                //					const long dof_I = slave_dofs[i];
// //                //
// //                //					for(int j = 0; j < master_dofs.size(); ++j) {
// //                //
// //                //                        const long dof_J = master_dofs[j];
// //                //
// //                //                        mat_buffer.add(dof_I, dof_J, elemmat(i, j));
// //                //					}
// //                //				}
// //
// //                return true;
// //
// //            } else {
// //
// //                return false;
// //            }
// //
// //        };
// //
// //        if(!SurfaceAssemble<Dimensions>(comm, master_slave, fun, settings, search_radius, tag_1, tag_2)) {
// //            return false;
// //        }
// //        //
// //        //
// //        //		// std::cout << mat_buffer << std::endl;
// //        //
// //        //		double volumes[2] = { local_element_matrices_sum,  total_intersection_volume };
// //        //
// //        //		comm.all_reduce(volumes, 2, moonolith::MPISum());
// //        //
// //        //		const processor_id_type master_proc_id  = master->mesh().processor_id();
// //        //
// //        //		const dof_id_type n_dofs_on_proc_master = master->dof_map().n_dofs_on_processor(master_proc_id);
// //        //
// //        //		const processor_id_type slave_proc_id   = slave->mesh().processor_id();
// //        //
// //        //		const dof_id_type n_dofs_on_proc_slave  = slave->dof_map().n_dofs_on_processor(slave_proc_id);
// //        //
// //        //		if(comm.is_root()) {
// //        //			std::cout << "sum(B): " << volumes[0] << ", vol(I): " << volumes[1] << std::endl;
// //        //		}
// //
// //        //
// //        //
// //        //		std::vector<moonolith::Integer>  ownershipRangesMaster(comm.size()+1);
// //        //		ownershipRangesMaster.allSet(0);
// //        //
// //        //
// //        //		std::vector<moonolith::Integer>  ownershipRangesSlave(comm.size()+1);
// //        //		ownershipRangesSlave.allSet(0);
// //        //
// //        //
// //        //		ownershipRangesMaster[comm.rank()+1]+= static_cast<unsigned int>(n_dofs_on_proc_master);
// //        //
// //        //		ownershipRangesSlave[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_slave);
// //        //
// //        //		comm.all_reduce(&ownershipRangesMaster[0], ownershipRangesMaster.size(), moonolith::MPIMax());
// //        //
// //        //		comm.all_reduce(&ownershipRangesSlave[0],  ownershipRangesSlave.size(),  moonolith::MPIMax());
// //        //
// //        //		moonolith::Redistribute< moonolith::SparseMatrix<double> > redist(comm.get_mpi_comm());
// //        //		redist.apply(ownershipRangesSlave, mat_buffer, moonolith::AddAssign<double>());
// //        //
// //        //		assert(ownershipRangesSlave.empty() == ownershipRangesMaster.empty() ||
// ownershipRangesMaster.empty());
// //        //
// //        //		moonolith::root_describe("petsc assembly begin", comm, std::cout);
// //        //
// //        //		SizeType  mMaxRowEntries = mat_buffer.local_max_entries_x_col();
// //        //		comm.all_reduce(&mMaxRowEntries, 1, moonolith::MPIMax());
// //        //
// //        //		const SizeType local_range_slave_range  = ownershipRangesSlave [comm.rank()+1] -
// ownershipRangesSlave [comm.rank()];
// //        //		const SizeType local_range_master_range = ownershipRangesMaster[comm.rank()+1] -
// ownershipRangesMaster[comm.rank()];
// //        //
// //        //		B = utopia::local_sparse(local_range_slave_range, local_range_master_range, mMaxRowEntries);
// //        //
// //        //		{
// //        //			utopia::Write<utopia::USparseMatrix> write(B);
// //        //			for (auto it = mat_buffer.iter(); it; ++it) {
// //        //				B.set(it.row(), it.col(), *it);
// //        //
// //        //			}
// //        //		}
// //        //
// //        //
// //        //		moonolith::root_describe("petsc assembly end", comm, std::cout);
// //        //
// //        //
// //        return true;
// //    }
// //    //
// //    //
// //    //
// //    //
// //
// //
// //
// //
// //    bool ParMortarAssembler::SurfaceAssemble(USparseMatrix &B, const libMesh::Real search_radius, const int tag_1,
// const int tag_2)
// //    {
// //
// //        moonolith::SearchSettings settings;
// //
// //
// //        moonolith::Communicator comm = libmesh_comm_.get();
// //
// //        if(master_slave_->mesh().mesh_dimension() == 2) {
// //            //std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
// //            return utopia::SurfaceAssemble<2>(comm, master_slave_, B, settings, search_radius, tag_1, tag_2);
// //        }
// //
// //
// //        if(master_slave_->mesh().mesh_dimension() == 3) {
// //            return utopia::SurfaceAssemble<3>(comm, master_slave_, B, settings, search_radius, tag_1, tag_2);
// //        }
// //
// //        assert(false && "Dimension not supported!");
// //        return false;
// //    }
// //    //
// //    //
// //    //
// //    //	bool ParMortarAssembler::Transfer(USparseMatrix &B, USparseMatrix &T)
// //    //	{
// //    //		UVector               diag_elem;
// //    //		PetscInt               mG;
// //    //		PetscInt               nG;
// //    //		PetscInt               mL;
// //    //		PetscInt               nL;
// //    //
// //    //		diag_elem = 1./sum(B,1);
// //    //		T = diag(diag_elem)*B;
// //    //
// //    //		return true;
// //    //	}
// //    //
// //}
// //
// //
// //

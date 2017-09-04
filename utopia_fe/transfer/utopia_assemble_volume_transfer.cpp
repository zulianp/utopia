#include "utopia_assemble_volume_transfer.hpp"

#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/serial_mesh.h"

#include "Box.hpp"
#include "utopia_fe_core.hpp"
#include "MortarAssemble.hpp"
#include "utopia_BoxAdapter.hpp"
#include "utopia_VElementAdapter.hpp"
#include "utopia_VTree.hpp"
#include "utopia_ElementDofMap.hpp"
#include "utopia_FESpacesAdapter.hpp"

#include "moonolith_profiler.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_tree.hpp"
#include "moonolith_n_tree_mutator_factory.hpp"
#include "moonolith_n_tree_with_span_mutator_factory.hpp"
#include "moonolith_n_tree_with_tags_mutator_factory.hpp"
#include "moonolith_sparse_matrix.hpp"
#include "par_moonolith.hpp"

#include <cmath>
#include <queue>


using namespace libMesh;

namespace utopia {
    
    template<typename T>
    inline void Print(const std::vector<T> &v, std::ostream &os)
    {
        for(auto i : v) {
            os << i << " ";
        }
        
        os << "\n";
    }

    
    template<class Iterator>
    static void write_space(const Iterator &begin, const Iterator &end,MeshBase &space,
                            const std::vector<ElementDofMap> &dof_map,/* const std::vector<ElementDofMap> &variable_number,*/
                            const std::vector<ElementDofMap> &variable_order,const int role, moonolith::OutputStream &os)
    {
        const int dim 		  = space.mesh_dimension();
        const long n_elements = std::distance(begin, end);
        
        std::set<long> nodeIds;
        std::map<long, long> mapping;
   
        for(Iterator it = begin; it != end; ++it) {
            
            const Elem *elem = space.elem(*it);
            
            for(dof_id_type j = 0; j != elem->n_nodes(); ++j) {
                
                nodeIds.insert(elem->node(j));
            }
            
        }
        
        long n_nodes = nodeIds.size();
        
        // Estimate for allocation
        os.request_space( (n_elements * 8 + n_nodes * dim) * (sizeof(double) + sizeof(long)) );
        
        //WRITE 1
        os << dim << role;
        
        
        int index = 0;
        for (auto nodeId : nodeIds) {
            mapping[nodeId] = index++;
        }
        
        //WRITE 2
        os << n_nodes;
        
        //WRITE 6
        os << n_elements;
        
        for(auto node_id : nodeIds){
            
            const Point &p = space.node(node_id);
            
            for(int i = 0; i < dim; ++i) {
                
                //WRITE 3
                os << p(i);
                
            }
            
            //std::cout<<"write_point"<<p<<std::endl;
            
        }
        
        std::vector<dof_id_type> indices_vector;
        
        CHECK_STREAM_WRITE_BEGIN("elements", os);
        
        for(Iterator it = begin; it != end; ++it) {
            
            const int k = *it;
            
            const Elem *elem = space.elem(*it);
            
            const int e_n_nodes = elem->n_nodes();
            
            const int type = elem->type();
            
            //WRITE 7
            os << type << e_n_nodes;
            
            
            
            for (int i = 0; i != e_n_nodes; ++i) {
                
                auto it = mapping.find(elem->node(i));
                
                assert(it != mapping.end());
                
                int index = it->second;
                
                //WRITE 8
                os << index;
                
            }
            
            //WRITE 9
            assert(!dof_map.at(elem->id()).empty());
            
            os << dof_map.at(elem->id());
            
            
            
        }
        
        CHECK_STREAM_WRITE_END("elements", os);
        
        //WRITE 10
        //        os << variable_number.at(0);
        
        //WRITE 11
        os << variable_order.at(0);
        
        
        
    }
    
    
    
    
    template<class Iterator>
    static void write_element_selection(const Iterator &begin, const Iterator &end, const FESpacesAdapter &spaces, moonolith::OutputStream &os)
    {
        
        if(spaces.spaces().empty()){
            assert(false);
            return;
        }
        
        auto m = spaces.spaces()[0];
        std::shared_ptr<MeshBase> s = nullptr;
        
        if(spaces.spaces().size()>1) {
            s=spaces.spaces()[1];
        }
        
        std::vector<long> master_selection;
        std::vector<long> slave_selection;
        
        bool met_slave_selection = false;
        
       
        for(Iterator it = begin; it != end; ++it) {
            int index =*it;
            
            if(m && index >= m->n_elem()) {
                index -= m->n_elem();
                slave_selection.push_back(index);
            }
            
            else if(!m) {
                met_slave_selection = true;
                slave_selection.push_back(index);
            }
            
            else {
                assert(!met_slave_selection);
                assert(index < m->n_elem());
                master_selection.push_back(index);
            }
        }
        
        
        
        const bool has_master = !master_selection.empty();
        const bool has_slave  = !slave_selection.empty();
        
        os << has_master << has_slave;
        
        
        if(has_master) {
            write_space(master_selection.begin(), master_selection.end(), *m, spaces.dof_map(0),
                        /*spaces.variable_number(0),*/ spaces.variable_order(0), 0, os);
        }
        
        if(has_slave) {
            write_space(slave_selection.begin(), slave_selection.end(), *s, spaces.dof_map(1),
                        /*spaces.variable_number(1),*/ spaces.variable_order(1), 1, os);
        }
        
    }

    static void read_space(moonolith::InputStream &is, std::shared_ptr<MeshBase> & space,
                           std::vector<ElementDofMap> &dof_map, /*std::vector<ElementDofMap> &variable_number,*/
                           std::vector<ElementDofMap> &variable_order, const libMesh::Parallel::Communicator &comm)
    {
     
        using namespace std;
        
        
        
        //READ 1
        int dim, role;
        is >> dim >> role;
        
        //READ 2
        long n_nodes;
        is >> n_nodes;
        
        //READ 6
        long n_elements;
        is >> n_elements;
        
        auto mesh_ptr = std::make_shared<SerialMesh>(comm, dim);
        
        mesh_ptr->reserve_nodes(n_nodes);
        
        for (long iii = 0; iii != n_nodes; ++iii) {
            
            Point p;
            
            for(int j = 0; j < dim; ++j) {
                //READ 3
                is >> p(j);
            }
            
            mesh_ptr->add_point(p);        
        }
        
        dof_map.resize(n_elements);
        
        
        CHECK_STREAM_READ_BEGIN("elements", is);
        
        for(long i = 0; i !=n_elements; ++i) {
            
            //READ 7
            
            int type, e_n_nodes;
            
            is >> type >> e_n_nodes;
            
            auto elem = Elem::build(ElemType(type)).release();
            
            
            int index;
            
            for (int ii = 0; ii != e_n_nodes; ++ii) {
                
                //READ 8
                is >> index;
                
                elem->set_node(ii) = & mesh_ptr->node(index);
                
            }
            
            //READ 9
            is >> dof_map.at(i);
            
            mesh_ptr->add_elem(elem);
            
            libmesh_assert(elem);
    
        }
        
        CHECK_STREAM_READ_END("elements", is);
        
        //READ 10
        //  variable_number.resize(1);
        //  is >> variable_number.at(0);
        
        //READ 11
        variable_order.resize(1);
        is >> variable_order.at(0);
        
        
        //!!!! dummy parameters
    
        space = mesh_ptr;
        
    }
    
    static void read_spaces(moonolith::InputStream &is, FESpacesAdapter &spaces, const libMesh::Parallel::Communicator &comm_master, const libMesh::Parallel::Communicator &comm_slave)
    {

        
        bool has_master, has_slave;
        is >> has_master >> has_slave;
        
        spaces.spaces().resize(2);
        
        
        if(has_master) {
            read_space(is, spaces.spaces()[0], spaces.dof_map(0),
                       /* spaces.variable_number(0),*/spaces.variable_order(0), comm_master);
            spaces.set_must_destroy_attached(0,true);
        } else {
            spaces.spaces()[0] = nullptr;
            
        }
        
        if(has_slave) {
            read_space(is, spaces.spaces()[1], spaces.dof_map(1),
                       /*spaces.variable_number(1),*/spaces.variable_order(1),comm_slave);
            spaces.set_must_destroy_attached(1,true);
        } else {
            spaces.spaces()[1] = nullptr;
            
        }
        
        
        
    }
    
    template<int Dimensions, class Fun>
    static bool Assemble(moonolith::Communicator &comm,
                         const std::shared_ptr<MeshBase> &master,
                         const std::shared_ptr<MeshBase> &slave,
                         const std::shared_ptr<DofMap> &dof_master,
                         const std::shared_ptr<DofMap> &dof_slave,
                         const unsigned int &_from_var_num,
                         const unsigned int &_to_var_num,
                         Fun process_fun,
                         const moonolith::SearchSettings &settings, 
                         bool use_biorth_, 
                         int n_var)
    {
        
        
        using namespace moonolith;
        
        typedef VTree<Dimensions> NTreeT;
        typedef typename NTreeT::DataContainer DataContainer;
        typedef typename NTreeT::DataType Adapter;
        
        const long maxNElements = settings.max_elements;
        const long maxDepth = settings.max_depth;
        
        
        const auto &master_mesh = master;
        const auto &slave_mesh  = slave;
        const int n_elements_master = master_mesh->n_elem();
        const int n_elements_slave  = slave_mesh->n_elem();
        const int n_elements 		= n_elements_master + n_elements_slave;
        
        
        const Parallel::Communicator &libmesh_comm_master = master_mesh->comm();
        const Parallel::Communicator &libmesh_comm_slave = slave_mesh->comm();
        
        
        auto predicate = std::make_shared<MasterAndSlave>();
        predicate->add(0, 1);
        
        MOONOLITH_EVENT_BEGIN("create_adapters");

        std::shared_ptr<NTreeT> tree = NTreeT::New(predicate, maxNElements, maxDepth);
        std::cout<<"n_elements"<<n_elements<<std::endl;
        tree->reserve(n_elements);
        
        
        auto local_spaces = std::make_shared<FESpacesAdapter>(master, slave, dof_master, dof_slave, _from_var_num, _to_var_num);
        
        int offset = 0;
        int space_num = 0;
        
        for(auto s : local_spaces->spaces()) {
            if(s) {
                
                bool first = true;
                for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it) {
                    auto elem=*it;
                    Adapter a(*s, elem->id(), offset+elem->id(), space_num);
                    assert(!local_spaces->dof_map(space_num)[elem->id()].empty());
                    a.set_dof_map(&local_spaces->dof_map(space_num)[elem->id()].global);
                    tree->insert(a);
                }
                
                offset += s->n_elem(); //s->mesh().n_active_local_elem();//(*s->mesh().active_local_elements_end())->id();
                
            }
            
            
            ++space_num;
            
            
        }
        
        tree->root()->bound().static_bound().enlarge(1e-8);
        
   
        MOONOLITH_EVENT_END("create_adapters");
        

        std::map<long, std::shared_ptr<FESpacesAdapter> > spaces;
        std::map<long, std::vector<std::shared_ptr<FESpacesAdapter> > > migrated_spaces;
        
        
        auto read = [&spaces, &migrated_spaces, comm, &libmesh_comm_master, &libmesh_comm_slave ]
        (
         const long ownerrank,
         const long senderrank,
         bool is_forwarding, DataContainer &data,
         InputStream &in
         ) {
            
        
            CHECK_STREAM_READ_BEGIN("vol_proj", in);
            
            std::shared_ptr<FESpacesAdapter> proc_space = std::make_shared<FESpacesAdapter>(comm);
            
            read_spaces(in, *proc_space, libmesh_comm_master, libmesh_comm_slave);
            
            if (!is_forwarding) {
                assert(!spaces[ownerrank]);
                spaces[ownerrank] = proc_space;
            } else {
                migrated_spaces[ownerrank].push_back(proc_space);
            }
            
            data.reserve(data.size() + 3000);

            
            int space_num = 0;
            long offset = 0;
            for(auto s : proc_space->spaces()) {
                if(s) {
                    for (int i=0; i<s->n_elem(); i++) {
                        data.push_back(Adapter(*s, i, offset + i, space_num) );
                        assert(!proc_space->dof_map(space_num)[i].empty());
                        data.back().set_dof_map(&proc_space->dof_map(space_num)[i].global);
                    }
                    
                    offset += s->n_elem();
                    
                }
                
                ++space_num;
                
            }
            
            
           
            CHECK_STREAM_READ_END("vol_proj", in);
            
            
            
        };
        
        
        auto write = [&local_spaces, &spaces, &comm]
        (
         const long ownerrank, const long recvrank,
         const std::vector<long>::const_iterator &begin,
         const std::vector<long>::const_iterator &end,
         const DataContainer &data,
         OutputStream &out) {
            
            CHECK_STREAM_WRITE_BEGIN("vol_proj", out);
            
         
            if (ownerrank == comm.rank()) {
                
                write_element_selection(begin, end, *local_spaces, out);
                
                
            } else {
                
                auto it = spaces.find(ownerrank);
                assert(it != spaces.end());
                std::shared_ptr<FESpacesAdapter> spaceptr = it->second;
                assert(std::distance(begin, end) > 0);
                write_element_selection(begin, end, *spaceptr, out);
                
            }

            
            CHECK_STREAM_WRITE_END("vol_proj", out);
            
        };
        
        
        long n_false_positives = 0, n_intersections = 0;
    
        
        auto fun = [&n_false_positives, &n_intersections, &process_fun](
                                                                        
                                                                        Adapter &master, Adapter &slave) -> bool {
            
            bool ok = process_fun(master, slave);
            
            if(ok) {
                n_intersections++;
                
                return true;
            } else {
                
                n_false_positives++;
                return false;
            }
            return true;
            
        };
         
        moonolith::search_and_compute(comm, tree, predicate, read, write, fun, settings);
        
        
        
        long n_total_candidates = n_intersections + n_false_positives;
        
        long n_collection[3] = {n_intersections, n_total_candidates, n_false_positives};
        
        
        comm.all_reduce(n_collection, 3, moonolith::MPISum());
        
        
        if (comm.is_root()) {
            std::cout << "n_intersections: " << n_collection[0]
            << ", n_total_candidates: " 	 << n_collection[1]
            << ", n_false_positives: " 	     << n_collection[2] << std::endl;
        }
        
        return true;
    }
    
    static void assemble_biorth_weights_from_space(const std::shared_ptr<MeshBase> &mesh,
                                                   const std::shared_ptr<DofMap> &dof_map,
                                                   const int var_num,
                                                   libMesh::DenseMatrix<libMesh::Real> &weights)
    {
        const int dim = mesh->mesh_dimension();
        std::unique_ptr<libMesh::FEBase> biorth_elem =
        libMesh::FEBase::build(dim,
                               dof_map->variable_type(var_num));
        
        auto &el = **mesh->active_local_elements_begin();
        
        const int order = order_for_l2_integral(dim,
                                                el, dof_map->variable(var_num).type().order,
                                                el, dof_map->variable(var_num).type().order);
        
        libMesh::QGauss qg(dim, libMesh::Order(order));
        biorth_elem->attach_quadrature_rule(&qg);
        biorth_elem->reinit(&el);
        mortar_assemble_weights(*biorth_elem, weights);
    }

    static void scale_polyhedron(const double scaling, Polyhedron &poly)
    {
       const int n_values = poly.n_nodes * poly.n_dims;
        for(int i = 0; i < n_values; ++i) {
            poly.points[i] *= scaling;
        }
    }
    
    template<int Dimensions>
    bool Assemble(
                  moonolith::Communicator &comm,
                  const std::shared_ptr<MeshBase> &master,
                  const std::shared_ptr<MeshBase> &slave,
                  const std::shared_ptr<DofMap> &dof_master,
                  const std::shared_ptr<DofMap> &dof_slave,
                  const unsigned int &_from_var_num,
                  const unsigned int &_to_var_num,
                  DSMatrixd &B,
                  const moonolith::SearchSettings &settings,bool  use_biorth_, int n_var)
    {
        
        const int var_num_slave = _to_var_num;
        
        std::shared_ptr<FESpacesAdapter> local_fun_spaces = std::make_shared<FESpacesAdapter>(master, slave, dof_master, dof_slave,_from_var_num,_to_var_num);
        
        libMesh::DenseMatrix<libMesh::Real> src_pts;
        libMesh::DenseMatrix<libMesh::Real> dest_pts;
        libMesh::DenseMatrix<libMesh::Real> intersection2;
        Polyhedron src_poly, dest_poly;
        Polyhedron  intersection3,temp_poly;
        Intersector isector;
        
        std::shared_ptr<MeshBase> master_space = master;
        std::shared_ptr<MeshBase> slave_space  = slave;

        libMesh::DenseMatrix<libMesh::Real> elemmat;
        libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;
        
        std::shared_ptr<Transform> src_trans;
        std::shared_ptr<Transform> dest_trans;
        
        
        int skip_zeros = 1;
        
        
        libMesh::Real total_intersection_volume = 0.0;
        libMesh::Real local_element_matrices_sum = 0.0;
        
        
        
        moonolith::SparseMatrix<double> mat_buffer(comm);
        mat_buffer.set_size(dof_slave->n_dofs(), dof_master->n_dofs());
        
        bool intersected = false;
        
        double element_setup_time = 0.0;
        double intersection_time = 0.0;
        double assembly_time     = 0.0;
        
        utopia::Chrono c;
        
        auto fun = [&](const VElementAdapter<Dimensions> &master,
                       const VElementAdapter<Dimensions> &slave) -> bool {
            
            c.start();
            
            libMesh::DenseMatrix<libMesh::Real> biorth_weights;
            
            if(use_biorth_) {
                assemble_biorth_weights_from_space(slave_space,
                                                   dof_slave,
                                                   var_num_slave,
                                                   biorth_weights);
            }
            
            long n_intersections = 0;
            
            bool pair_intersected = false;
            
            const auto &src  = master.space();
            
            const auto &dest = slave.space();
            
            const auto &src_mesh  = src;
            
            const auto &dest_mesh = dest;
            
            const int src_index  = master.element();
            
            const int dest_index = slave.element();
            
            auto &src_el  = *src_mesh.elem(src_index);
            
            auto &dest_el = *dest_mesh.elem(dest_index);
            
            const int dim = src_mesh.mesh_dimension();
            
            
            std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
            
            master_fe = libMesh::FEBase::build(src_mesh.mesh_dimension(),  dof_master->variable_type(0));
            slave_fe  = libMesh::FEBase::build(dest_mesh.mesh_dimension(), dof_slave->variable_type(0));
            
            QMortar composite_ir(dim);
            QMortar src_ir(dim);
            QMortar dest_ir(dim);
            
            
            const int order = order_for_l2_integral(dim, src_el, dof_master->variable(0).type().order , dest_el,dof_slave->variable(0).type().order);
            
            c.stop();
            element_setup_time += c.get_seconds();
            c.start();
            
            if(dim == 2)  {
                make_polygon(src_el,   src_pts);
                make_polygon(dest_el, dest_pts);
                
                if(intersect_2D(src_pts, dest_pts, intersection2)) {
                    total_intersection_volume += fabs(isector.polygon_area_2(intersection2.m(), &intersection2.get_values()[0]));
                    
                    const libMesh::Real weight=isector.polygon_area_2(dest_pts.m(), &dest_pts.get_values()[0]);
                    
                    make_composite_quadrature_2D(intersection2, weight, order, composite_ir);
                    pair_intersected = true;
                    //bool affine_map = true;
                    src_trans  = std::make_shared<AffineTransform2>(src_el);
                    dest_trans = std::make_shared<AffineTransform2>(dest_el);
                    pair_intersected = true;
                }
            }
            else if(dim == 3) {

                // const libMesh::Real weight = isector.p_mesh_volume_3(dest_poly); //(volume)
                //src_poly
                make_polyhedron(src_el,  src_poly);
                make_polyhedron(dest_el, dest_poly);

                //scale_polyhedron(src_poly, 1./weight); //(volume)
                //scale_polyhedron(dest_poly, 1./weight); //(volume)
                
                
                if(intersect_3D(src_poly, dest_poly, intersection3)) {
                    //scale_polyhedron(intersection3, weight); //(volume)

                    total_intersection_volume += isector.p_mesh_volume_3(intersection3);
                    
                    const libMesh::Real weight = isector.p_mesh_volume_3(dest_poly); //comment this out (volume)
                    
                    make_composite_quadrature_3D(intersection3, weight, order, composite_ir);
                    src_trans  = std::make_shared<AffineTransform3>(src_el);
                    dest_trans = std::make_shared<AffineTransform3>(dest_el);
                    pair_intersected = true;
                }
                
            } else {
                assert(false);
                return false;
            }
            
            c.stop();
            intersection_time += c.get_seconds();
            c.start();
            
            const auto &master_dofs = master.dof_map();
            const auto &slave_dofs  = slave.dof_map();
            
            if(pair_intersected) {
                
                
                transform_to_reference(*src_trans,  src_el.type(),  composite_ir,  src_ir);
                transform_to_reference(*dest_trans, dest_el.type(), composite_ir,  dest_ir);
                

                assert(!master_dofs.empty());
                assert(!slave_dofs.empty());

                master_fe->attach_quadrature_rule(&src_ir);

                const std::vector<std::vector<Real>> & phi_master  = master_fe->get_phi();

                master_fe->reinit(&src_el);
                
                slave_fe->attach_quadrature_rule(&dest_ir);

                const std::vector<std::vector<Real>> & phi_slave = slave_fe->get_phi();
           
                const std::vector<Real> & JxW_slave = slave_fe->get_JxW();

                slave_fe->reinit(&dest_el);
                
                elemmat.zero();
                
                if(use_biorth_) {
                    mortar_assemble_weighted_biorth(*master_fe, *slave_fe, biorth_weights, elemmat);
                } else {
                    mortar_assemble(*master_fe, *slave_fe, elemmat);
                }
            
                auto partial_sum = std::accumulate(elemmat.get_values().begin(), elemmat.get_values().end(), libMesh::Real(0.0));
 
                local_element_matrices_sum += partial_sum;
                
                intersected = true;
                
                ++n_intersections;
                

                assert(slave_dofs.size() == elemmat.m());
                assert(master_dofs.size() == elemmat.n());
                
                         
                for(int i = 0; i < slave_dofs.size(); ++i) {
                    
                    const long dof_I = slave_dofs[i];
                    
                    for(int j = 0; j < master_dofs.size(); ++j) {
                        
                        const long dof_J = master_dofs[j];
                        
                        mat_buffer.add(dof_I, dof_J, elemmat(i, j));
                    }
                }
                
                return true;
                
            } else {
                
                return false;
            }
            
        };

        if(!Assemble<Dimensions>(comm, master, slave, dof_master, dof_slave, _from_var_num, _to_var_num, fun, settings, use_biorth_,n_var)) {
            std::cout << "no intersections" <<std::endl;
            return false;
        }
                
        double volumes[2] = { local_element_matrices_sum,  total_intersection_volume };
        
        comm.all_reduce(volumes, 2, moonolith::MPISum());
        
        const processor_id_type master_proc_id  = master->processor_id();
        
        const dof_id_type n_dofs_on_proc_master = dof_master->n_local_dofs();
        
        const processor_id_type slave_proc_id   = slave->processor_id();
        
        const dof_id_type n_dofs_on_proc_slave  =dof_slave->n_local_dofs();
        
        const int n_dofs_on_proc_print  = dof_slave->n_local_dofs();
        
        if(comm.is_root()) {
            std::cout << "sum(B): " << volumes[0] << ", vol(I): " << volumes[1] << std::endl;
        }
        
        std::vector<moonolith::Integer> ownershipRangesMaster(comm.size()+1, 0);        
        std::vector<moonolith::Integer> ownershipRangesSlave(comm.size()+1, 0);

        ownershipRangesMaster[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_master);        
        ownershipRangesSlave[comm.rank()+1]  += static_cast<unsigned int>(n_dofs_on_proc_slave);
        
        comm.all_reduce(&ownershipRangesMaster[0], ownershipRangesMaster.size(), moonolith::MPISum());
        comm.all_reduce(&ownershipRangesSlave[0],  ownershipRangesSlave.size(),  moonolith::MPISum());
        
        std::partial_sum(ownershipRangesMaster.begin(), ownershipRangesMaster.end(), ownershipRangesMaster.begin());
        std::partial_sum(ownershipRangesSlave.begin(), ownershipRangesSlave.end(), ownershipRangesSlave.begin());
        
        
        int dim = master->mesh_dimension();
        
        moonolith::Redistribute< moonolith::SparseMatrix<double> > redist(comm.get_mpi_comm());
        
        redist.apply(ownershipRangesSlave, mat_buffer, moonolith::AddAssign<double>());
        
        assert(ownershipRangesSlave.empty() == ownershipRangesMaster.empty() || ownershipRangesMaster.empty());
        
        moonolith::root_describe("petsc assembly begin", comm, std::cout);
        
        SizeType  mMaxRowEntries = mat_buffer.local_max_entries_x_col();
        
        comm.all_reduce(&mMaxRowEntries, 1, moonolith::MPIMax());
        
        const SizeType local_range_slave_range  = ownershipRangesSlave [comm.rank()+1] - ownershipRangesSlave [comm.rank()];
        const SizeType local_range_master_range = ownershipRangesMaster[comm.rank()+1] - ownershipRangesMaster[comm.rank()];
        
        DSMatrixd B_x = utopia::local_sparse(local_range_slave_range, local_range_master_range, mMaxRowEntries);
        
        {
            utopia::Write<utopia::DSMatrixd> write(B_x);
            for (auto it = mat_buffer.iter(); it; ++it) {
                B_x.set(it.row(), it.col(), *it);
                
            }
        }
        
        
        auto s_B_x = local_size(B_x);
        
        B = local_sparse(s_B_x.get(0), s_B_x.get(1), n_var * mMaxRowEntries);
        
        utopia::Write<DSMatrixd> w_B(B);
        utopia::each_read(B_x, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
            for(utopia::SizeType d = 0; d < n_var; ++d) {
                B.set(i+d, j+d, value);
            }
        });
        
        moonolith::root_describe("petsc assembly end", comm, std::cout);
        return true;
        
    }
    
    
    
    bool assemble_volume_transfer(moonolith::Communicator &comm,
                              const std::shared_ptr<MeshBase> &master,
                              const std::shared_ptr<MeshBase> &slave,
                              const std::shared_ptr<DofMap> &dof_master,
                              const std::shared_ptr<DofMap> &dof_slave,
                              const unsigned int & _from_var_num,
                              const unsigned int & _to_var_num,
                              bool  use_biorth_, int n_var, DSMatrixd &B)
    {
        moonolith::SearchSettings settings;
        
        if(master->mesh_dimension() == 2) {
            return utopia::Assemble<2>(comm, master, slave, dof_master, dof_slave, _from_var_num,  _to_var_num, B, settings,use_biorth_, n_var);
        }
        
        
        if(master->mesh_dimension() == 3) {
            return utopia::Assemble<3>(comm, master, slave, dof_master, dof_slave, _from_var_num,  _to_var_num, B, settings,use_biorth_, n_var);
        }
        
        assert(false && "Dimension not supported!");
        return false;
    }
}




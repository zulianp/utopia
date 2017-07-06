#include "utopia_assemble_volume_transfer_r.hpp"

#include "cutlibpp.hpp"
#include "cutlibpp_Base.hpp"
#include "cutlibpp_Tree.hpp"
#include "cutlibpp_NTreeMutatorFactory.hpp"
#include "cutlibpp_NTreeWithSpanMutatorFactory.hpp"
#include "cutlibpp_NTreeWithTagsMutatorFactory.hpp"
#include "cutlibpp_API.hpp"

#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/serial_mesh.h"

#include "express_Profiler.hpp"
#include "express_Redistribute.hpp"

#include "Box.hpp"
#include "MapSparseMatrix.hpp"
#include "utopia_fe.hpp"
#include "MortarAssemble.hpp"
#include "libmesh/serial_mesh.h"
#include "utopia_BoxAdapter.hpp"
#include "utopia_VElementAdapter.hpp"
#include "utopia_VTree.hpp"
#include "utopia_ElementDofMap.hpp"
#include "utopia_FESpacesRAdapter.hpp"
#include <cmath>
#include <queue>



namespace utopia {
    using namespace libMesh;
    
    template<typename T>
    inline void Print(const std::vector<T> &v, std::ostream &os)
    {
        for(auto i : v) {
            os << i << " ";
        }
        
        os << "\n";
    }
    
    static std::ostream &logger()
    {
        return express::Express::Instance().logger().os();
    }
    
    
    template<class Iterator>
    static void write_space(const Iterator &begin, const Iterator &end,MeshBase &space,
                            const std::vector<ElementDofMap> &dof_map,/* const std::vector<ElementDofMap> &variable_number,*/
                            const std::vector<ElementDofMap> &dof_map_reverse,
                            const std::vector<ElementDofMap> &variable_order,const int role, cutk::OutputStream &os)
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
        os.requestSpace( (n_elements * 8 + n_nodes * dim) * (sizeof(double) + sizeof(long)) );
        
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
        
        //        std::cout<<"write_n_el = "<<n_elements<<std::endl;
        
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
            assert(!dof_map_reverse.at(elem->id()).empty());
            
            os << dof_map.at(elem->id());
            os << dof_map_reverse.at(elem->id());
            
            
            
        }
        
        CHECK_STREAM_WRITE_END("elements", os);
        
        //WRITE 10
        //        os << variable_number.at(0);
        
        //WRITE 11
        os << variable_order.at(0);
        
        
        
    }
    
    
    
    
    template<class Iterator>
    static void write_element_selection(const Iterator &begin, const Iterator &end, const FESpacesRAdapter &spaces, cutk::OutputStream &os)
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
            write_space(master_selection.begin(), master_selection.end(), *m, spaces.dof_map(0), spaces.dof_map_reverse(0),
                        /*spaces.variable_number(0),*/ spaces.variable_order(0), 0, os);
        }
        
        if(has_slave) {
            write_space(slave_selection.begin(), slave_selection.end(), *s, spaces.dof_map(1), spaces.dof_map_reverse(1),
                        /*spaces.variable_number(1),*/ spaces.variable_order(1), 1, os);
        }
        
        
    }
    
    static void read_space(cutk::InputStream &is, cutk::shared_ptr<MeshBase> & space,
                           std::vector<ElementDofMap> &dof_map, /*std::vector<ElementDofMap> &variable_number,*/
                           std::vector<ElementDofMap> &dof_map_reverse,
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
        
        dof_map_reverse.resize(n_elements);
        
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
            is >> dof_map_reverse.at(i);
            
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
    
    static void read_spaces(cutk::InputStream &is, FESpacesRAdapter &spaces, const libMesh::Parallel::Communicator &comm_master, const libMesh::Parallel::Communicator &comm_slave)
    {
        
        bool has_master, has_slave;
        is >> has_master >> has_slave;
        
        spaces.spaces().resize(2);
        
        
        if(has_master) {
            read_space(is, spaces.spaces()[0], spaces.dof_map(0),spaces.dof_map_reverse(0),
                       /* spaces.variable_number(0),*/spaces.variable_order(0), comm_master);
            spaces.set_must_destroy_attached(0,true);
        } else {
            spaces.spaces()[0] = nullptr;
            
        }
        
        if(has_slave) {
            read_space(is, spaces.spaces()[1], spaces.dof_map(1), spaces.dof_map_reverse(1),
                       /*spaces.variable_number(1),*/spaces.variable_order(1),comm_slave);
            spaces.set_must_destroy_attached(1,true);
        } else {
            spaces.spaces()[1] = nullptr;
            
        }
    }
    
    template<int Dimensions, class Fun>
    static bool Assemble(express::Communicator &comm,
                         const std::shared_ptr<MeshBase> &master,
                         const std::shared_ptr<MeshBase> &slave,
                         const std::shared_ptr<DofMap> &dof_master,
                         const std::shared_ptr<DofMap> &dof_slave,
                         const std::shared_ptr<DofMap> &dof_reverse_master,
                         const std::shared_ptr<DofMap> &dof_reverse_slave,
                         const unsigned int &_from_var_num,
                         const unsigned int &_to_var_num,
                         const unsigned int &_from_var_num_r,
                         const unsigned int &_to_var_num_r,
                         Fun process_fun,
                         const cutk::Settings &settings, bool use_biorth_, int n_var, int n_var_r)
    {
        
        
        using namespace cutlibpp;
        using namespace express;
        using namespace cutk;
        
        typedef VTree<Dimensions> NTreeT;
        typedef typename NTreeT::DataContainer DataContainer;
        typedef typename NTreeT::DataType Adapter;
        
        long maxNElements = 100;
        long maxDepth = 5;
        
        
        if (!settings.get("max_depth").isNull()) {
            maxDepth = settings.get("max_depth").toInt();
            std::cout<<"max_depth  = "<< maxDepth  <<std::endl;
        }
        
        const auto &master_mesh = master;
        const auto &slave_mesh  = slave;
        const int n_elements_master = master_mesh->n_elem();
        const int n_elements_slave  = slave_mesh->n_elem();
        const int n_elements 		= n_elements_master + n_elements_slave;
        
        
        const Parallel::Communicator &libmesh_comm_master = master_mesh->comm();
        const Parallel::Communicator &libmesh_comm_slave = slave_mesh->comm();
        
        
        auto predicate = make_shared<MasterAndSlave>();
        predicate->add(0, 1);
        
        EXPRESS_EVENT_BEGIN("create_adapters");
  
        cutk::shared_ptr<NTreeT> tree = NTreeT::New(predicate, maxNElements, maxDepth);
        std::cout<<"n_elements"<<n_elements<<std::endl;
        tree->reserve(n_elements);
        
        
        std::shared_ptr<FESpacesRAdapter> local_spaces = make_shared<FESpacesRAdapter>(master, slave, dof_master, dof_slave, dof_reverse_master, dof_reverse_slave, _from_var_num, _to_var_num, _from_var_num_r, _to_var_num_r);
        int offset = 0;
        int space_num = 0;
        
        for(auto s : local_spaces->spaces()) {
            if(s) {
                
                bool first = true;
                for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it) {
                    auto elem=*it;
                    Adapter a(*s, elem->id(), offset+elem->id(), space_num);
                    assert(!local_spaces->dof_map(space_num)[elem->id()].empty());
                    assert(!local_spaces->dof_map_reverse(space_num)[elem->id()].empty());
                    a.set_dof_map(&local_spaces->dof_map(space_num)[elem->id()].global);
                    a.set_dof_map_reverse(&local_spaces->dof_map_reverse(space_num)[elem->id()].global);
                    tree->insert(a);
                }
                
                offset += s->n_elem(); //s->mesh().n_active_local_elem();//(*s->mesh().active_local_elements_end())->id();
                
            }
            
            
            ++space_num;
            
            
        }
        
        std::cout<<" --------------------------- tree->memory().nData()=" << tree->memory().nData()<<std::endl; ;
        
        tree->getRoot()->getBound().staticBound().enlarge(1e-8);
        

        EXPRESS_EVENT_END("create_adapters");
        
        //Just to have an indexed-storage
        std::map<long, cutk::shared_ptr<FESpacesRAdapter> > spaces;
        std::map<long, std::vector<cutk::shared_ptr<FESpacesRAdapter> > > migrated_spaces;
        
        
        auto read = [&spaces, &migrated_spaces, comm, &libmesh_comm_master, &libmesh_comm_slave ]
        (
         const long ownerrank,
         const long senderrank,
         bool is_forwarding, DataContainer &data,
         InputStream &in
         ) {
            
            CHECK_STREAM_READ_BEGIN("vol_proj", in);
            
            cutk::shared_ptr<FESpacesRAdapter> proc_space = cutk::make_shared<FESpacesRAdapter>(comm);
            
            read_spaces(in, *proc_space, libmesh_comm_master, libmesh_comm_slave);
            
            if (!is_forwarding) {
                assert(!spaces[ownerrank]);
                spaces[ownerrank] = proc_space;
            } else {
                migrated_spaces[ownerrank].push_back(proc_space);
            }
            
            data.reserve(data.size() + 3000);
            //std::cout<<proc_space->dof_map(0)[0].global<<std::endl;
            
            int space_num = 0;
            long offset = 0;
            for(auto s : proc_space->spaces()) {
                if(s) {
                    for (int i=0; i<s->n_elem(); i++) {
                        data.push_back(Adapter(*s, i, offset + i, space_num) );
                        assert(!proc_space->dof_map(space_num)[i].empty());
                        assert(!proc_space->dof_map_reverse(space_num)[i].empty());
                        data.back().set_dof_map(&proc_space->dof_map(space_num)[i].global);
                        data.back().set_dof_map_reverse(&proc_space->dof_map_reverse(space_num)[i].global);
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
                cutk::shared_ptr<FESpacesRAdapter> spaceptr = it->second;
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
        
        
        
        cutk::Settings custom_settings = settings;
        custom_settings.set("disable_redistribution", cutk::Boolean(true));
        custom_settings.set("verbosity_level", cutk::Integer(1));
        
        
        cutlibpp::search_and_compute(comm, tree, predicate, read, write, fun, custom_settings);
        
        
        
        long n_total_candidates = n_intersections + n_false_positives;
        
        long n_collection[3] = {n_intersections, n_total_candidates, n_false_positives};
        
        
        comm.allReduce(n_collection, 3, express::MPISum());
        
        
        if (comm.isRoot()) {
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
    
    template<int Dimensions>
    bool Assemble(
                  express::Communicator &comm,
                  const std::shared_ptr<MeshBase> &master,
                  const std::shared_ptr<MeshBase> &slave,
                  const std::shared_ptr<DofMap> &dof_master,
                  const std::shared_ptr<DofMap> &dof_slave,
                  const std::shared_ptr<DofMap> &dof_reverse_master,
                  const std::shared_ptr<DofMap> &dof_reverse_slave,
                  const unsigned int &_from_var_num,
                  const unsigned int &_to_var_num,
                  const unsigned int &_from_var_num_r,
                  const unsigned int &_to_var_num_r,
                  DSMatrixd &B, DSMatrixd &B_reverse, //bbecsek
                  const cutk::Settings &settings,bool  use_biorth_, int n_var, int n_var_r)
    {
        
        const int var_num_slave = _to_var_num;
        
        std::shared_ptr<FESpacesRAdapter> local_fun_spaces = cutk::make_shared<FESpacesRAdapter>(master, slave, dof_master, dof_slave, dof_reverse_master, dof_reverse_slave, _from_var_num, _to_var_num, _from_var_num_r,
        _to_var_num_r);
        
        libMesh::DenseMatrix<libMesh::Real> master_pts;
        libMesh::DenseMatrix<libMesh::Real> slave_pts;
        libMesh::DenseMatrix<libMesh::Real> intersection2;
        Polyhedron master_poly, slave_poly;
        Polyhedron  intersection3,temp_poly;
        Intersector isector;
        
        std::shared_ptr<MeshBase> master_space = master;
        std::shared_ptr<MeshBase> slave_space  = slave;
        
        
        // std::vector<libMesh::dof_id_type> master_dofs, slave_dofs;
        libMesh::DenseMatrix<libMesh::Real> elemmat;
        libMesh::DenseMatrix<libMesh::Real> elemmat_reverse;
        libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;
        
        std::shared_ptr<Transform> master_trans;
        std::shared_ptr<Transform> slave_trans;
        
        
        int skip_zeros = 1;
        
        
        libMesh::Real total_intersection_volume = 0.0;
        libMesh::Real local_element_matrices_sum = 0.0;
        libMesh::Real local_element_matrices_sum_reverse = 0.0;
        
        
        
        
        express::MapSparseMatrix<double> mat_buffer(dof_slave->n_dofs(), dof_master->n_dofs());

        express::MapSparseMatrix<double> mat_buffer_reverse(dof_reverse_master->n_dofs(), dof_reverse_slave->n_dofs());
  
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
            
            const auto &master_mesh  = src;
            
            const auto &slave_mesh = dest;
            
            const int master_index  = master.element();
            
            const int slave_index = slave.element();
            
            auto &master_el  = *master_mesh.elem(master_index);
            
            auto &slave_el = *slave_mesh.elem(slave_index);
            
            const int dim = master_mesh.mesh_dimension();
            
            
            std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
            
            master_fe = libMesh::FEBase::build(master_mesh.mesh_dimension(),  dof_master->variable_type(0));
            slave_fe  = libMesh::FEBase::build(slave_mesh.mesh_dimension(), dof_slave->variable_type(0));
            
            QMortar composite_ir(dim);
            QMortar master_ir(dim);
            QMortar slave_ir(dim);
            
            
            const int order = order_for_l2_integral(dim, master_el, dof_master->variable(0).type().order , slave_el,dof_slave->variable(0).type().order);
            libMesh::Real weight_reverse = 0;
            
            c.stop();
            element_setup_time += c.get_seconds();
            c.start();
            
            if(dim == 2)  {
                make_polygon(master_el,   master_pts);
                make_polygon(slave_el, slave_pts);
                
                if(intersect_2D(master_pts, slave_pts, intersection2)) {
                    total_intersection_volume += fabs(isector.polygon_area_2(intersection2.m(), &intersection2.get_values()[0]));
                    
                    const libMesh::Real weight = isector.polygon_area_2(slave_pts.m(), &slave_pts.get_values()[0]);
                    weight_reverse = isector.polygon_area_2(master_pts.m(), &master_pts.get_values()[0])/weight;
                    
                    make_composite_quadrature_2D(intersection2, weight, order, composite_ir);
                    pair_intersected = true;
                    
                    master_trans  = std::make_shared<Transform2>(master_el);
                    slave_trans = std::make_shared<Transform2>(slave_el);
                    pair_intersected = true;
                }
            }
            else if(dim == 3) {
                make_polyhedron(master_el,  master_poly);
                make_polyhedron(slave_el, slave_poly);
                
                
                if(intersect_3D(master_poly, slave_poly, intersection3)) {
                    
                    total_intersection_volume += isector.p_mesh_volume_3(intersection3);
                    
                    const libMesh::Real weight = isector.p_mesh_volume_3(slave_poly);
                    weight_reverse = isector.p_mesh_volume_3(master_poly)/weight;
                    
                    make_composite_quadrature_3D(intersection3, weight, order, composite_ir);
                    master_trans  = std::make_shared<Transform3>(master_el);
                    slave_trans = std::make_shared<Transform3>(slave_el);
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
            
            const auto &master_dofs_reverse = master.dof_map_reverse();
            const auto &slave_dofs_reverse  = slave.dof_map_reverse();
            
            if(pair_intersected) {
                
                
                transform_to_reference(*master_trans,  master_el.type(),  composite_ir,  master_ir);
                transform_to_reference(*slave_trans, slave_el.type(), composite_ir,  slave_ir);
                
                //important for correct scaling of the quadrature weights
                
                for(int i = 0; i < master_ir.n_points(); ++i) {
                    master_ir.get_weights()[i] /= weight_reverse;
                }
              
                assert(!master_dofs.empty());
                assert(!slave_dofs.empty());
            
                master_fe->attach_quadrature_rule(&master_ir);
                master_fe->get_phi();
                master_fe->get_JxW();
                master_fe->reinit(&master_el);
                
                slave_fe->attach_quadrature_rule(&slave_ir);
                slave_fe->get_phi();
                slave_fe->get_JxW();
                slave_fe->reinit(&slave_el);
                
                elemmat.zero();
                elemmat_reverse.zero();
                
                
                if(use_biorth_) {
                    mortar_assemble_weighted_biorth(*master_fe, *slave_fe, biorth_weights, elemmat);

                    mortar_assemble_weighted_biorth(*slave_fe, *master_fe, biorth_weights, elemmat);
                    
                } else {
                    mortar_assemble(*master_fe, *slave_fe, elemmat);

                    mortar_assemble(*slave_fe, *master_fe, elemmat_reverse);
                }
                
                
                local_element_matrices_sum +=         std::accumulate(elemmat.get_values().begin(), elemmat.get_values().end(), libMesh::Real(0.0));
                local_element_matrices_sum_reverse += std::accumulate(elemmat_reverse.get_values().begin(), elemmat_reverse.get_values().end(), libMesh::Real(0.0));;
                
                intersected = true;
                
                ++n_intersections;
                           
                assert(slave_dofs.size() == elemmat.m());
                assert(master_dofs.size() == elemmat.n());

                // std::cout<<"master_dofs_reverse.size() ==>"<< master_dofs_reverse.size()<<std::endl;

                // std::cout<<"elemmat_reverse.m() ==>"<< elemmat_reverse.m()<<std::endl;
                
                assert(master_dofs_reverse.size() == elemmat_reverse.m());
                assert(slave_dofs_reverse.size() == elemmat_reverse.n());

                for(int i = 0; i < slave_dofs.size(); ++i) {
                    
                    const long dof_I = slave_dofs[i];
                    
                    for(int j = 0; j < master_dofs.size(); ++j) {
                        
                        const long dof_J = master_dofs[j];
                        
                        mat_buffer.add(dof_I, dof_J, elemmat(i, j));
                    }
                }
                
//                 bbecsek
                for(int i = 0; i < master_dofs_reverse.size(); ++i) {
                    
                    const long dof_I_r = master_dofs_reverse[i];
            
                    for(int j = 0; j < slave_dofs_reverse.size(); ++j) {
                        
                        const long dof_J_r = slave_dofs_reverse[j];
                        
                        mat_buffer_reverse.add(dof_I_r, dof_J_r, elemmat_reverse(i, j));
                    }
                }

                return true;
                
            } else {
                
                return false;
            }
            
        };
        
        

        if(!Assemble<Dimensions>(comm, master, slave, dof_master, dof_slave, dof_reverse_master, dof_reverse_slave, _from_var_num, _to_var_num, _from_var_num_r, _to_var_num_r, fun, settings, use_biorth_, n_var, n_var_r)) {

            return false;
        }
  
        
        
        double volumes[3] = { local_element_matrices_sum,  total_intersection_volume, local_element_matrices_sum_reverse };
        
        comm.allReduce(volumes, 3, express::MPISum());
        
        const processor_id_type master_proc_id  = master->processor_id();
        
        const dof_id_type n_dofs_on_proc_master = dof_master->n_local_dofs();
        
        const processor_id_type slave_proc_id   = slave->processor_id();
        
        const dof_id_type n_dofs_on_proc_slave  =dof_slave->n_local_dofs();
        
        const int n_dofs_on_proc_print  = dof_slave->n_local_dofs();
        
        
        const dof_id_type n_dofs_on_proc_master_r  = dof_reverse_master->n_local_dofs();
        const dof_id_type n_dofs_on_proc_slave_r   = dof_reverse_slave->n_local_dofs();
        
        std::cout<<" dof_slave_r->n_local_dofs() " <<  dof_reverse_master->n_local_dofs() <<std::endl;
        
        std::cout<<" dof_master_r->n_local_dofs() " << dof_reverse_slave->n_local_dofs() <<std::endl;
        
        
        if(comm.isRoot()) {
            std::cout << "sum(B): " << volumes[0] << ", vol(I): " << volumes[1] << std::endl;
            std::cout << "sum(B*): " << volumes[2] << std::endl;
        }
        
        express::Array<express::SizeType>  ownershipRangesMaster(comm.size()+1);
        ownershipRangesMaster.allSet(0);
        
        express::Array<express::SizeType>  ownershipRangesSlave(comm.size()+1);
        ownershipRangesSlave.allSet(0);
        
        
        ownershipRangesMaster[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_master);
        
        ownershipRangesSlave[comm.rank()+1]  += static_cast<unsigned int>(n_dofs_on_proc_slave);
        
        
        
        comm.allReduce(&ownershipRangesMaster[0], ownershipRangesMaster.size(), express::MPISum());
        
        comm.allReduce(&ownershipRangesSlave[0],  ownershipRangesSlave.size(),  express::MPISum());
        
        
        std::partial_sum(ownershipRangesMaster.begin(), ownershipRangesMaster.end(), ownershipRangesMaster.begin());
        std::partial_sum(ownershipRangesSlave.begin(), ownershipRangesSlave.end(), ownershipRangesSlave.begin());
        
        
        // bbecsek:
        express::Array<express::SizeType>  ownershipRangesMaster_r(comm.size()+1);
        ownershipRangesMaster_r.allSet(0);
        
        express::Array<express::SizeType>  ownershipRangesSlave_r(comm.size()+1);
        ownershipRangesSlave_r.allSet(0);
        
        ownershipRangesMaster_r[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_master_r);
        
        ownershipRangesSlave_r[comm.rank()+1]  += static_cast<unsigned int>(n_dofs_on_proc_slave_r);
        
        comm.allReduce(&ownershipRangesMaster_r[0], ownershipRangesMaster_r.size(), express::MPISum());
        
        comm.allReduce(&ownershipRangesSlave_r[0],  ownershipRangesSlave_r.size(),  express::MPISum());
        
        std::partial_sum(ownershipRangesMaster_r.begin(), ownershipRangesMaster_r.end(), ownershipRangesMaster_r.begin());
        std::partial_sum(ownershipRangesSlave_r.begin(), ownershipRangesSlave_r.end(), ownershipRangesSlave_r.begin());
        //
        
        
        //        if(comm.isRoot()) {
        //            std::cout <<ownershipRangesMaster << std::endl;
        //            std::cout<<"prova"<<n_dofs_on_proc_print<<std::endl;
        //
        //        }
        
        
        
        int dim = master->mesh_dimension();
        
        express::Redistribute< express::MapSparseMatrix<double> > redist(comm.getMPIComm());
        
        // bbecsek
        express::Redistribute< express::MapSparseMatrix<double> > redist_reverse(comm.getMPIComm());
        
        redist.apply(ownershipRangesSlave, mat_buffer, express::AddAssign<double>());
        
        // bbecsek
        redist_reverse.apply(ownershipRangesMaster_r, mat_buffer_reverse, express::AddAssign<double>());
        
        assert(ownershipRangesSlave.empty() == ownershipRangesMaster.empty() || ownershipRangesMaster.empty());
        
        assert(ownershipRangesMaster_r.empty() == ownershipRangesSlave_r.empty() || ownershipRangesSlave_r.empty());
        
        express::RootDescribe("petsc assembly begin", comm, std::cout);
        
        SizeType  mMaxRowEntries = mat_buffer.maxEntriesXCol();
        
        // bbecsek
        SizeType  mMaxRowEntries_reverse = mat_buffer_reverse.maxEntriesXCol();
        
        comm.allReduce(&mMaxRowEntries, 1, express::MPIMax());
        comm.allReduce(&mMaxRowEntries_reverse, 1, express::MPIMax());
        
        const SizeType local_range_slave  = ownershipRangesSlave [comm.rank()+1] - ownershipRangesSlave [comm.rank()];
        const SizeType local_range_master = ownershipRangesMaster[comm.rank()+1] - ownershipRangesMaster[comm.rank()];
        // bbecsek
        
        express::RootDescribe("petsc assembly begin 1", comm, std::cout);
        const SizeType local_range_slave_r  = ownershipRangesSlave_r [comm.rank()+1] - ownershipRangesSlave_r [comm.rank()];
        const SizeType local_range_master_r = ownershipRangesMaster_r[comm.rank()+1] - ownershipRangesMaster_r[comm.rank()];
        
        express::RootDescribe("petsc assembly begin 2", comm, std::cout);
        
        DSMatrixd B_x = utopia::local_sparse(local_range_slave, local_range_master, mMaxRowEntries);
        // bbecsek
        DSMatrixd B_x_reverse = utopia::local_sparse(local_range_master_r, local_range_slave_r, mMaxRowEntries_reverse);
//        
        {
            utopia::Write<utopia::DSMatrixd> write(B_x);
            for (auto it = mat_buffer.iter(); it; ++it) {
                B_x.set(it.row(), it.col(), *it);
                
            }
        }
        
        
        utopia::write("B_x.m",B_x);
        // bbecsek
        
        express::RootDescribe("petsc assembly begin 3", comm, std::cout);
        {
            utopia::Write<utopia::DSMatrixd> write_reverse(B_x_reverse);
            for (auto it = mat_buffer_reverse.iter(); it; ++it) {
                B_x_reverse.set(it.row(), it.col(), *it);
                
            }
        }
        
        //utopia::write("B_x_reverse.m",B_x_reverse);
        
        auto s_B_x = local_size(B_x);
        // bbecsek
        auto s_B_x_reverse = local_size(B_x_reverse);
        
        B = local_sparse(s_B_x.get(0), s_B_x.get(1), n_var * mMaxRowEntries);
        // bbecsek
        B_reverse = local_sparse(s_B_x_reverse.get(0), s_B_x_reverse.get(1), n_var_r * mMaxRowEntries_reverse);
        
        
        std::cout<< "modify the matrix  B"<<std::endl;
        utopia::Write<DSMatrixd> w_B(B);
        utopia::each_read(B_x, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
            for(utopia::SizeType d = 0; d < n_var; ++d) {
                B.set(i+d, j+d, value);
            }
        });
        
        //bbecsek: can we move this to the first iteration?
        std::cout<<"n_var_r"<<n_var_r<<std::endl;
        std::cout<< "modify the matrix  B_reverse"<<std::endl;
        utopia::Write<DSMatrixd> w_B_reverse(B_reverse);
        utopia::each_read(B_x_reverse, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
            for(utopia::SizeType d = 0; d < n_var_r ; ++d) {
                B_reverse.set(i+d, j+d, value);
            }
        });
        
        
        //disp(B_reverse.size());
        
        
        
        
        express::RootDescribe("petsc assembly end", comm, std::cout);
        
        return true;
        
    }
    
    
    
    
    bool assemble_volume_transfer_r(express::Communicator &comm,
                                    const std::shared_ptr<MeshBase> &master,
                                    const std::shared_ptr<MeshBase> &slave,
                                    const std::shared_ptr<DofMap> &dof_master,
                                    const std::shared_ptr<DofMap> &dof_slave,
                                    const std::shared_ptr<DofMap> &dof_reverse_master,
                                    const std::shared_ptr<DofMap> &dof_reverse_slave,
                                    const unsigned int & _from_var_num,
                                    const unsigned int & _to_var_num,
                                    const unsigned int & _from_var_num_r,
                                    const unsigned int & _to_var_num_r,
                                    bool  use_biorth_,
                                    int n_var,
                                    int n_var_r,
                                    DSMatrixd &B,
                                    DSMatrixd &B_reverse)
    {
        cutk::Settings settings;
        
        if(master->mesh_dimension() == 2) {
            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
            return utopia::Assemble<2>(comm,
                                       master, slave,
                                       dof_master, dof_slave,
                                       dof_reverse_master, dof_reverse_slave,
                                       _from_var_num,  _to_var_num,
                                       _from_var_num_r, _to_var_num_r,
                                       B, B_reverse,
                                       settings,use_biorth_, n_var, n_var_r);
        }
        
        
        if(master->mesh_dimension() == 3) {
            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
            return utopia::Assemble<3>(comm,
                                       master, slave,
                                       dof_master, dof_slave,
                                       dof_reverse_master, dof_reverse_slave,
                                       _from_var_num,  _to_var_num,
                                       _from_var_num_r, _to_var_num_r,
                                       B, B_reverse,
                                       settings,use_biorth_, n_var, n_var_r);
        }
        
        assert(false && "Dimension not supported!");
        return false;
    }

}




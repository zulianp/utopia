#include "utopia_assemble_volume_transfer_r.hpp"

#include "Box.hpp"
#include "utopia_fe_core.hpp"
#include "MortarAssemble.hpp"
#include "utopia_BoxAdapter.hpp"
#include "utopia_VElementAdapter.hpp"
#include "utopia_VTree.hpp"
#include "utopia_ElementDofMap.hpp"
#include "utopia_FESpacesRAdapter.hpp"

#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/serial_mesh.h"

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
#include <algorithm>
#include <unordered_set>

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
    
    template<int Dimensions, class Fun>
    static bool Assemble(moonolith::Communicator &comm,
                         const std::shared_ptr<MeshBase> &master,
                         const std::shared_ptr<MeshBase> &slave,
                         const std::shared_ptr<DofMap> &dof_master,
                         const std::shared_ptr<DofMap> &dof_slave,
                         const std::shared_ptr<DofMap> &dof_reverse_master,
                         const std::shared_ptr<DofMap> &dof_reverse_slave,
                         const unsigned int &from_var_num,
                         const unsigned int &to_var_num,
                         const unsigned int &from_var_num_r,
                         const unsigned int &to_var_num_r,
                         Fun process_fun,
                         const moonolith::SearchSettings &settings, 
                         bool use_biorth_, 
                         int n_var, 
                         int n_var_r)
    {
        using namespace moonolith;
        
        typedef VTree<Dimensions> NTreeT;
        typedef typename NTreeT::DataContainer DataContainer;
        typedef typename NTreeT::DataType Adapter;
        
        const long maxNElements = settings.max_elements;
        const long maxDepth = settings.max_depth;
                
        const auto &master_mesh = master;
        const auto &slave_mesh  = slave;
        const int n_elements_master = master_mesh->n_active_local_elem();
        const int n_elements_slave  = slave_mesh->n_active_local_elem();
        const int n_elements 		= n_elements_master + n_elements_slave;
        
        
        const Parallel::Communicator &libmesh_comm_master = master_mesh->comm();
        const Parallel::Communicator &libmesh_comm_slave = slave_mesh->comm();
        
        auto predicate = std::make_shared<MasterAndSlave>();
        predicate->add(0, 1);
        
        MOONOLITH_EVENT_BEGIN("create_adapters");
  
        auto tree = NTreeT::New(predicate, maxNElements, maxDepth);
        tree->reserve(n_elements);
        
        auto local_spaces = std::make_shared<FESpacesRAdapter>(master, slave, dof_master, dof_slave, dof_reverse_master, dof_reverse_slave, from_var_num, to_var_num, from_var_num_r, to_var_num_r);
        int offset = 0;
        int space_num = 0;
        
        for(auto s : local_spaces->spaces()) {
            if(s) {
                bool first = true;
                libMesh::dof_id_type local_element_id = 0;
                for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it, ++local_element_id) {
                    auto elem = *it;
                    Adapter a(*s, elem->id(), offset + local_element_id, space_num);
                    assert(!local_spaces->dof_map(space_num)[local_element_id].empty());
                    assert(!local_spaces->dof_map_reverse(space_num)[local_element_id].empty());
                    a.set_dof_map(&local_spaces->dof_map(space_num)[local_element_id].global);
                    a.set_dof_map_reverse(&local_spaces->dof_map_reverse(space_num)[local_element_id].global);
                    tree->insert(a);
                }
                
                offset += s->n_active_local_elem();
            }
            
            ++space_num;
        }
                
        tree->root()->bound().static_bound().enlarge(1e-8);
        
        MOONOLITH_EVENT_END("create_adapters");
        
        //Just to have an indexed-storage
        std::map<long, std::shared_ptr<FESpacesRAdapter> > spaces;
        std::map<long, std::vector<std::shared_ptr<FESpacesRAdapter> > > migrated_spaces;

        auto read = [&spaces, &migrated_spaces, comm, &libmesh_comm_master, &libmesh_comm_slave]
        (
         const long ownerrank,
         const long senderrank,
         bool is_forwarding, DataContainer &data,
         InputStream &in
         ) {
            
            CHECK_STREAM_READ_BEGIN("vol_proj", in);
            
            std::shared_ptr<FESpacesRAdapter> proc_space = std::make_shared<FESpacesRAdapter>(comm);
            
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
                std::shared_ptr<FESpacesRAdapter> spaceptr = it->second;
                assert(std::distance(begin, end) > 0);
                write_element_selection(begin, end, *spaceptr, out);
            }

            CHECK_STREAM_WRITE_END("vol_proj", out);
            
        };
        
        
        long n_false_positives = 0, n_intersections = 0;
        
        auto fun = [&n_false_positives, &n_intersections, &process_fun](Adapter &master, Adapter &slave) -> bool {
            
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
    
    template<int Dimensions>
    bool Assemble(
                  moonolith::Communicator &comm,
                  const std::shared_ptr<MeshBase> &master,
                  const std::shared_ptr<MeshBase> &slave,
                  const std::shared_ptr<DofMap> &dof_master,
                  const std::shared_ptr<DofMap> &dof_slave,
                  const std::shared_ptr<DofMap> &dof_reverse_master,
                  const std::shared_ptr<DofMap> &dof_reverse_slave,
                  const unsigned int &from_var_num,
                  const unsigned int &to_var_num,
                  const unsigned int &from_var_num_r,
                  const unsigned int &to_var_num_r,
                  DSMatrixd &B, 
                  DSMatrixd &B_reverse, //bbecsek
                  const moonolith::SearchSettings &settings,
                  std::unordered_set<utopia::SizeType> &particip_slave_dofs,
                  bool use_biorth_, 
                  int n_var, 
                  int n_var_r)
    {
        
        const int var_num_slave = to_var_num;
        
        auto local_fun_spaces = std::make_shared<FESpacesRAdapter>(master, slave, dof_master, dof_slave, dof_reverse_master, dof_reverse_slave, from_var_num, to_var_num, from_var_num_r,
        to_var_num_r);
        
        libMesh::DenseMatrix<libMesh::Real> master_pts;
        libMesh::DenseMatrix<libMesh::Real> slave_pts;
        libMesh::DenseMatrix<libMesh::Real> intersection2;
        Polyhedron master_poly, slave_poly;
        Polyhedron  intersection3,temp_poly;
        Intersector isector;
        
        std::shared_ptr<MeshBase> master_space = master;
        std::shared_ptr<MeshBase> slave_space  = slave;
        
        libMesh::DenseMatrix<libMesh::Real> elemmat;
        libMesh::DenseMatrix<libMesh::Real> elemmat_reverse;
        libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;
        
        std::shared_ptr<Transform> master_trans;
        std::shared_ptr<Transform> slave_trans;
        
        int skip_zeros = 1;
        
        
        libMesh::Real total_intersection_volume = 0.0;
        libMesh::Real local_element_matrices_sum = 0.0;
        libMesh::Real local_element_matrices_sum_reverse = 0.0;
        
        moonolith::SparseMatrix<double> mat_buffer(comm);
        mat_buffer.set_size(dof_slave->n_dofs(), dof_master->n_dofs());
       
        moonolith::SparseMatrix<double> mat_buffer_reverse(comm);
        mat_buffer_reverse.set_size(dof_reverse_master->n_dofs(), dof_reverse_slave->n_dofs());
  
        bool intersected = false;
        
        double element_setup_time = 0.0;
        double intersection_time = 0.0;
        double assembly_time     = 0.0;
        
        //utopia::Chrono c;
        
        auto fun = [&](const VElementAdapter<Dimensions> &master,
                       const VElementAdapter<Dimensions> &slave) -> bool {
            
            //c.start();
            
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
            
            //c.stop();
            //element_setup_time += c.get_seconds();
            //c.start();
            
            if(dim == 2)  {
                make_polygon(master_el,   master_pts);
                make_polygon(slave_el, slave_pts);
                
                //for(unsigned int node = 0; node < slave_el.n_nodes(); ++node)
                //  std::cout << slave_el.node_id(node) << std::endl; 

                if(intersect_2D(master_pts, slave_pts, intersection2)) {
                    total_intersection_volume += fabs(isector.polygon_area_2(intersection2.m(), &intersection2.get_values()[0]));
                    
                    const libMesh::Real weight = isector.polygon_area_2(slave_pts.m(), &slave_pts.get_values()[0]);
                    weight_reverse = isector.polygon_area_2(master_pts.m(), &master_pts.get_values()[0])/weight;
                    
                    make_composite_quadrature_2D(intersection2, weight, order, composite_ir);
                    pair_intersected = true;
                    
//                    bool affine_transf=true;
                    
                    master_trans  = std::make_shared<AffineTransform2>(master_el);
                    slave_trans = std::make_shared<AffineTransform2>(slave_el);
                    pair_intersected = true;
                    //debug
             //       for(unsigned int node = 0; node < slave_el.n_nodes(); ++node){
             //         auto &nodep = slave_el.node_ref(node);
             //         auto x = nodep(0); 
             //         auto y = nodep(1); 
             //         auto z = nodep(2);
             //         FILE *coords = fopen("candidate_coords.txt","a");
             //         fprintf(coords, "%8.4E,%8.4E,%8.4E\n", x, y, z);
             //         fclose(coords);  
             //         if (slave_el.node_id(node) == 7046)
             //           std::cout << "WTF!!!! This node should not be intersecting!!!! FTW\n";
             //       }
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
                    master_trans  = std::make_shared<AffineTransform3>(master_el);
                    slave_trans = std::make_shared<AffineTransform3>(slave_el);
                    pair_intersected = true;
                }
                
            } else {
                assert(false);
                return false;
            }
            
//            c.stop();
//            intersection_time += c.get_seconds();
//            c.start();
            
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
        
        

        if(!Assemble<Dimensions>(comm, master, slave, dof_master, dof_slave, dof_reverse_master, dof_reverse_slave, from_var_num, to_var_num, from_var_num_r, to_var_num_r, fun, settings, use_biorth_, n_var, n_var_r)) {

            return false;
        }
  
        
        
        double volumes[3] = { local_element_matrices_sum,  total_intersection_volume, local_element_matrices_sum_reverse };
        
        comm.all_reduce(volumes, 3, moonolith::MPISum());
        
        const processor_id_type master_proc_id  = master->processor_id();
        
        const dof_id_type n_dofs_on_proc_master = dof_master->n_local_dofs();
        
        const processor_id_type slave_proc_id   = slave->processor_id();
        
        const dof_id_type n_dofs_on_proc_slave  =dof_slave->n_local_dofs();
        
        const int n_dofs_on_proc_print  = dof_slave->n_local_dofs();
        
        
        const dof_id_type n_dofs_on_proc_master_r  = dof_reverse_master->n_local_dofs();
        const dof_id_type n_dofs_on_proc_slave_r   = dof_reverse_slave->n_local_dofs();
        
        std::cout<<" dof_slave_r->n_local_dofs() " <<  dof_reverse_master->n_local_dofs() <<std::endl;
        
        std::cout<<" dof_master_r->n_local_dofs() " << dof_reverse_slave->n_local_dofs() <<std::endl;
        
        
        if(comm.is_root()) {
            std::cout << "sum(B): " << volumes[0] << ", vol(I): " << volumes[1] << std::endl;
            std::cout << "sum(B*): " << volumes[2] << std::endl;
        }
        
        std::vector<moonolith::Integer>  ownershipRangesMaster(comm.size()+1, 0);        
        std::vector<moonolith::Integer>  ownershipRangesSlave(comm.size()+1, 0);        
        
        ownershipRangesMaster[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_master);
        ownershipRangesSlave[comm.rank()+1]  += static_cast<unsigned int>(n_dofs_on_proc_slave);
        
        comm.all_reduce(&ownershipRangesMaster[0], ownershipRangesMaster.size(), moonolith::MPISum());
        comm.all_reduce(&ownershipRangesSlave[0],  ownershipRangesSlave.size(),  moonolith::MPISum());
        
        
        std::partial_sum(ownershipRangesMaster.begin(), ownershipRangesMaster.end(), ownershipRangesMaster.begin());
        std::partial_sum(ownershipRangesSlave.begin(), ownershipRangesSlave.end(), ownershipRangesSlave.begin());
        
        
        // bbecsek:
        std::vector<moonolith::Integer>  ownershipRangesMaster_r(comm.size()+1, 0);
        std::vector<moonolith::Integer>  ownershipRangesSlave_r(comm.size()+1, 0);
        
        ownershipRangesMaster_r[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_master_r);
        ownershipRangesSlave_r[comm.rank()+1]  += static_cast<unsigned int>(n_dofs_on_proc_slave_r);
        
        comm.all_reduce(&ownershipRangesMaster_r[0], ownershipRangesMaster_r.size(), moonolith::MPISum());
        comm.all_reduce(&ownershipRangesSlave_r[0],  ownershipRangesSlave_r.size(),  moonolith::MPISum());
        
        std::partial_sum(ownershipRangesMaster_r.begin(), ownershipRangesMaster_r.end(), ownershipRangesMaster_r.begin());
        std::partial_sum(ownershipRangesSlave_r.begin(), ownershipRangesSlave_r.end(), ownershipRangesSlave_r.begin());
 
        int dim = master->mesh_dimension();
        
        moonolith::Redistribute< moonolith::SparseMatrix<double> > redist(comm.get_mpi_comm());
        
        // bbecsek
        moonolith::Redistribute< moonolith::SparseMatrix<double> > redist_reverse(comm.get_mpi_comm());
        
        redist.apply(ownershipRangesSlave, mat_buffer, moonolith::AddAssign<double>());
        
        // bbecsek
        redist_reverse.apply(ownershipRangesMaster_r, mat_buffer_reverse, moonolith::AddAssign<double>());
        
        assert(ownershipRangesSlave.empty() == ownershipRangesMaster.empty() || ownershipRangesMaster.empty());
        
        assert(ownershipRangesMaster_r.empty() == ownershipRangesSlave_r.empty() || ownershipRangesSlave_r.empty());
        
        moonolith::root_describe("petsc assembly begin", comm, std::cout);
        
        SizeType  mMaxRowEntries = mat_buffer.max_entries_x_col();
        
        // bbecsek
        SizeType  mMaxRowEntries_reverse = mat_buffer_reverse.max_entries_x_col();
        
        comm.all_reduce(&mMaxRowEntries, 1, moonolith::MPIMax());
        comm.all_reduce(&mMaxRowEntries_reverse, 1, moonolith::MPIMax());
        
        const SizeType local_range_slave  = ownershipRangesSlave [comm.rank()+1] - ownershipRangesSlave [comm.rank()];
        const SizeType local_range_master = ownershipRangesMaster[comm.rank()+1] - ownershipRangesMaster[comm.rank()];
        // bbecsek
        
        moonolith::root_describe("petsc assembly begin 1", comm, std::cout);
        const SizeType local_range_slave_r  = ownershipRangesSlave_r [comm.rank()+1] - ownershipRangesSlave_r [comm.rank()];
        const SizeType local_range_master_r = ownershipRangesMaster_r[comm.rank()+1] - ownershipRangesMaster_r[comm.rank()];
        
        moonolith::root_describe("petsc assembly begin 2", comm, std::cout);
        
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
        
        
        // utopia::write("B_x.m",B_x);
        // bbecsek
        
        moonolith::root_describe("petsc assembly begin 3", comm, std::cout);
        
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
        
        
        // std::cout<< "modify the matrix  B"<<std::endl;
        particip_slave_dofs.clear();
        utopia::Write<DSMatrixd> w_B(B);
        utopia::each_read(B_x, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
            for(utopia::SizeType d = 0; d < n_var; ++d) {
                B.set(i+d, j+d, value);
                particip_slave_dofs.emplace(i+d); // here we collect all dofs of the slave grid to be able to 
                                                  // restrict our mass matrix
            }
        });
        
        //bbecsek: can we move this to the first iteration?
        // std::cout<<"n_var_r"<<n_var_r<<std::endl;
        // std::cout<< "modify the matrix  B_reverse"<<std::endl;
        utopia::Write<DSMatrixd> w_B_reverse(B_reverse);
        utopia::each_read(B_x_reverse, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
            for(utopia::SizeType d = 0; d < n_var_r ; ++d) {
                B_reverse.set(i+d, j+d, value);
            }
        });
                
        //int world_rank;
        //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        //std::cout << "Participate size: " << particip_slave_dofs.size() << " on processor " << world_rank << std::endl;
 
        //disp(B_reverse.size());
        moonolith::root_describe("petsc assembly end", comm, std::cout);
        return true;
    }
    
    
    
    
    bool assemble_volume_transfer_r(moonolith::Communicator &comm,
                                    const std::shared_ptr<MeshBase> &master,
                                    const std::shared_ptr<MeshBase> &slave,
                                    const std::shared_ptr<DofMap> &dof_master,
                                    const std::shared_ptr<DofMap> &dof_slave,
                                    const std::shared_ptr<DofMap> &dof_reverse_master,
                                    const std::shared_ptr<DofMap> &dof_reverse_slave,
                                    const unsigned int & from_var_num,
                                    const unsigned int & to_var_num,
                                    const unsigned int & from_var_num_r,
                                    const unsigned int & to_var_num_r,
                                    std::unordered_set<utopia::SizeType> &particip_slave_dofs,
                                    bool  use_biorth_,
                                    int n_var,
                                    int n_var_r,
                                    DSMatrixd &B,
                                    DSMatrixd &B_reverse)
    {
        moonolith::SearchSettings settings;
        
        if(master->mesh_dimension() == 2) {
            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
            return utopia::Assemble<2>(comm,
                                       master, 
                                       slave,
                                       dof_master, 
                                       dof_slave,
                                       dof_reverse_master, 
                                       dof_reverse_slave,
                                       from_var_num,  
                                       to_var_num,
                                       from_var_num_r, 
                                       to_var_num_r,
                                       B, 
                                       B_reverse,
                                       settings,
                                       particip_slave_dofs,
                                       use_biorth_, 
                                       n_var, 
                                       n_var_r);
        }
        
        
        if(master->mesh_dimension() == 3) {
            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
            return utopia::Assemble<3>(comm,
                                       master, 
                                       slave,
                                       dof_master, 
                                       dof_slave,
                                       dof_reverse_master, 
                                       dof_reverse_slave,
                                       from_var_num,  
                                       to_var_num,
                                       from_var_num_r, 
                                       to_var_num_r,
                                       B, 
                                       B_reverse,
                                       settings,
                                       particip_slave_dofs,
                                       use_biorth_, 
                                       n_var, 
                                       n_var_r);
        }
        
        assert(false && "Dimension not supported!");
        return false;
    }

}




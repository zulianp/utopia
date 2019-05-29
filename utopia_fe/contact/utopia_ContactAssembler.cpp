#include "utopia_ContactAssembler.hpp"

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

#include "moonolith_affine_transform.hpp"
#include "moonolith_contact.hpp"
#include "moonolith_sparse_matrix.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_assign_functions.hpp"
#include "moonolith_elem_shape.hpp"
#include "moonolith_elem_triangle.hpp"
#include "moonolith_elem_quad.hpp"
#include "moonolith_elem_segment.hpp"

#include <vector>
#include <memory>

#include "libmesh/mesh_refinement.h"

namespace utopia {

    
    void ContactTensors::convert(const USparseMatrix &perm,
                                 const USparseMatrix &vector_perm,
                                 ContactTensors &out) const
    {
        out.B = vector_perm * B * transpose(vector_perm);
        out.D = vector_perm * D * transpose(vector_perm);
        
        if(!empty(Q)) {
            out.Q = vector_perm * Q * transpose(vector_perm);
            normalize_rows(out.Q);

            out.Q_inv = vector_perm * Q_inv * transpose(vector_perm);
            normalize_rows(out.Q_inv);
        }
        
        out.weighted_gap    = perm * weighted_gap;
        out.weighted_normal = vector_perm * weighted_normal;
        out.is_contact      = perm * is_contact;
        out.is_glue         = perm * is_glue;
    
        SizeType n_local_contact_nodes = 0;
        each_transform(out.is_contact, out.is_contact, [&n_local_contact_nodes](const SizeType i, const double val) -> double {
            if(val > 0) {
                ++n_local_contact_nodes;
                return 1.;
            } else {
                return 0.0;
            }
        });

        //remove filted out nodes (glue nodes are a subset of contact nodes)
        out.is_glue = e_mul(out.is_glue, out.is_contact);

        SizeType n_local_glue_nodes = 0;
        each_transform(out.is_glue, out.is_glue, [&n_local_glue_nodes](const SizeType i, const double val) -> double {
            if(val > 1e-8) {
                ++n_local_glue_nodes;
                return 1.;
            } else {
                return 0.0;
            }
        });

        double sum_normal   = sum(out.weighted_normal);
        double sum_normal_e = sum(weighted_normal);
        
        std::cout << "n_local_contact_nodes: " << n_local_contact_nodes << std::endl;
        std::cout << "n_local_glue_nodes:    " << n_local_glue_nodes << std::endl;
        std::cout << "sum(normal):           " << sum_normal << " == " << sum_normal_e << std::endl;
    }
    
    bool ContactTensors::check_op(const USparseMatrix &T)
    {
        UVector sum_T = sum(T, 1);
        
        bool ok = true;
        each_read(sum_T, [&ok](const SizeType i, const double val) {
            if(!approxeq(val, 0.0, 1e-1) && !approxeq(val, 1.0, 1e-1)) {
                std::cerr << i << "] " << val << "\n";
                ok = false;
            }
        });
        
        return ok;
    }
    
    void ContactTensors::finalize(const SizeType spatial_dim, const bool normalize)
    {
        zero_rows_to_identity(D, 1e-13);
        gap = local_zeros(local_size(weighted_gap));
        
        if(!empty(Q)) {
            UVector d_inv = diag(D);

            e_pseudo_inv(d_inv, d_inv, 1e-15);

            inv_mass_vector = d_inv;
            
            USparseMatrix D_tilde_inv = diag(d_inv);
            D_inv = Q * D_tilde_inv;
            T = D_inv * B;
            
            normalize_rows(T);
            
            assert(check_op(T));

            T += local_identity(local_size(T));
            
            gap    = Q * e_mul(d_inv, weighted_gap);
            normal = Q * e_mul(d_inv, weighted_normal);
            
        } else {

            inv_mass_vector = diag(D);
            e_pseudo_inv(inv_mass_vector, inv_mass_vector, 1e-15);

            gap    = e_mul(inv_mass_vector, weighted_gap);
            normal = e_mul(inv_mass_vector, weighted_normal);

            D_inv = diag(inv_mass_vector);
            T = D_inv * B;
            
            normalize_rows(T);
            assert(check_op(T));

            T += local_identity(local_size(T));
        }
        
        if(normalize) {
            auto r = range(normal);
            
            Read<UVector> ric(is_contact);
            ReadAndWrite<UVector> rw(normal);
            for(auto i = r.begin(); i < r.end(); i += spatial_dim) {
                bool is_node_in_contact = is_contact.get(i);
                
                if(is_node_in_contact) {
                    double len_n = 0.0;
                    for(SizeType d = 0; d < spatial_dim; ++d) {
                        auto temp = normal.get(i + d);
                        len_n += temp * temp;
                    }
                    
                    if(len_n > 0.0) {
                        len_n = std::sqrt(len_n);
                        
                        for(SizeType d = 0; d < spatial_dim; ++d) {
                            normal.set(i + d, normal.get(i + d)/len_n);
                        }
                    } else {
                        assert(false);
                    }
                } else {
                    for(SizeType d = 0; d < spatial_dim; ++d) {
                        normal.set(i + d, 0.0);
                    }
                }
            }
        }
        
        auto ls = local_size(D);
        orthogonal_trafo = local_sparse(ls.get(0), ls.get(1), spatial_dim);
        
        {
            auto r = range(normal);
            
            Read<UVector> ric(is_contact), rn(normal);
            Write<USparseMatrix> w(orthogonal_trafo);
            
            
            std::vector<double> nn(spatial_dim, 0.);
            OrthogonalTransformation ot;
            libMesh::DenseMatrix<double> local_ot;
            
            for(auto i = r.begin(); i < r.end(); i += spatial_dim) {
                bool is_node_in_contact = is_contact.get(i);
                
                if(is_node_in_contact) {
                    for(SizeType d = 0; d < spatial_dim; ++d) {
                        nn[d] = normal.get(i + d);
                    }
                    
                    ot.build(nn, local_ot);
                    
                    for(SizeType di = 0; di < spatial_dim; ++di) {
                        for(SizeType dj = 0; dj < spatial_dim; ++dj) {
                            orthogonal_trafo.set(i + di, i + dj, local_ot(di, dj));
                        }
                    }
                } else {
                    
                    for(SizeType di = 0; di < spatial_dim; ++di) {
                        orthogonal_trafo.set(i + di, i + di, 1.0);
                    }
                }
            }
        }


        {
            auto r = range(is_glue);
            ReadAndWrite<UVector> rw(is_glue);
            for(auto i = r.begin(); i < r.end(); i += spatial_dim) {
                auto val = is_glue.get(i);

                if(val > 0.) {
                    for(SizeType d = 1; d < spatial_dim; ++d) {
                        is_glue.set(i + d, val);
                    }
                }
            }
        }
        
        complete_transformation = T * orthogonal_trafo;
    }

    // void ContactTensors::determine_glued(const ContactParams &params)
    // {
    //     is_glue_node = local_zeros(local_size(is_contact));
    // }
    
    
    class ContactDataBuffers {
    public:
        
        inline MPI_Comm comm() const
        {
            return B.comm.get();
        }
        
        void finalize(const LibMeshFunctionSpaceAdapter &adapter,
                      const UVector &volumes,
                      const double alpha,
                      const bool use_biorth)
        {
            const SizeType n_local_elems = adapter.n_local_elems();
            const SizeType n_local_dofs  = adapter.n_local_dofs();
            const SizeType spatial_dim   = adapter.spatial_dim();


            
            area.finalize(n_local_elems);
            area.fill(element_wise.area);
            
            auto r = range(volumes);
            
            std::vector<bool> remove(adapter.n_local_dofs(), false);
            auto cr = adapter.permutation()->implementation().col_range();
            
            element_wise.is_contact = local_zeros(adapter.n_local_dofs());
            
            UVector elem_to_transform;
            
            const bool must_assemble_trafo = use_biorth && adapter.fe_type(0).order == 2;
            elem_to_transform = local_zeros(r.extent());
            
            {
                Read<UVector> rv(volumes), ra(element_wise.area);
                Write<UVector> wic(element_wise.is_contact), wett(elem_to_transform);
                
                for(auto i = r.begin(); i < r.end(); ++i) {
                    const auto a = element_wise.area.get(i);
                    
                    if(a > 0.0) {
                        auto ratio = a/volumes.get(i);
                        
                        if(!approxeq(ratio, 1.0, 1e-2)) {
                            remove[i - r.begin()] = true;
                            
                            const auto &dofs = adapter.element_dof_map()[i - r.begin()].global;
                            
                            for(auto d : dofs) {
                                remove[d - cr.begin()] = true;
                            }

                            std::cout << "=====================================\n";
                            std::cout << i << ") removed. " << ratio << " % of slave volume" << std::endl;
                            
                        } else {
                            // std::cout << "=====================================\n";
                            // std::cout << i << ") " << volumes.get(i) << " == " << element_wise.area.get(i) << std::endl;
                            
                            const auto &dofs = adapter.element_dof_map()[i - r.begin()].global;
                            
                            for(auto d : dofs) {
                                element_wise.is_contact.set(d, 1.0);
                            }
                            
                            if(must_assemble_trafo) {
                                elem_to_transform.set(i, 1.);
                            }
                        }
                    }
                }
            }
            
            B.finalize(n_local_dofs, n_local_dofs);
            D.finalize(n_local_dofs, n_local_dofs);
            gap.finalize(n_local_dofs);
            normal.finalize(n_local_dofs * spatial_dim);
            is_glue.finalize(n_local_dofs);
            
            B.fill(remove, element_wise.B_x);
            D.fill(remove, element_wise.D_x);
            gap.fill(remove, element_wise.weighted_gap);
            normal.fill(element_wise.weighted_normal);
            is_glue.fill(element_wise.is_glue);
            
            tensor_prod_with_identity(element_wise.B_x, spatial_dim, element_wise.B);
            tensor_prod_with_identity(element_wise.D_x, spatial_dim, element_wise.D);
            
            if(must_assemble_trafo) {
                DualBasis::build_global_trafo(
                                              adapter.mesh(),
                                              adapter.n_local_dofs(),
                                              adapter.element_dof_map(),
                                              elem_to_transform,
                                              alpha,
                                              element_wise.Q_x,
                                              false
                                              );

                
                tensor_prod_with_identity(element_wise.Q_x, spatial_dim, element_wise.Q);


                DualBasis::build_global_trafo(
                                              adapter.mesh(),
                                              adapter.n_local_dofs(),
                                              adapter.element_dof_map(),
                                              elem_to_transform,
                                              alpha,
                                              element_wise.Q_x,
                                              true);
                
                
                tensor_prod_with_identity(element_wise.Q_x, spatial_dim, element_wise.Q_inv);
            }
            
            double sum_B_x = sum(element_wise.B_x);
            double sum_D_x = sum(element_wise.D_x);
            double sum_normal_xyz = sum(element_wise.weighted_normal);
            
            assert(adapter.permutation());
            
            element_wise.convert(*adapter.permutation(),
                                 *adapter.vector_permutation(),
                                 *dof_wise
                                 );
            
            dof_wise->finalize(spatial_dim);
            
            
            static const double LARGE_VALUE = 1e6;
            
            Read<UVector> ric(dof_wise->is_contact);
            each_transform(dof_wise->gap, dof_wise->gap, [&](const SizeType i, const double value) -> double {
                const SizeType offset = i;
                
                if(dof_wise->is_contact.get(i) > 0) {
                    return value;
                } else {
                    return LARGE_VALUE;
                }
            });
        }
        
        MatrixInserter B, D;
        MatrixInserter gap, normal;
        MatrixInserter area;
        MatrixInserter is_glue;
        ContactTensors element_wise;
        std::shared_ptr<ContactTensors> dof_wise;
        
        ContactDataBuffers(MPI_Comm comm) : B(comm), D(comm), gap(comm), normal(comm), area(comm), is_glue(comm, false) {
            dof_wise = std::make_shared<ContactTensors>();
        }
    private:
        ContactDataBuffers(const ContactDataBuffers &other) 
        : B(other.comm()), D(other.comm()), gap(other.comm()), normal(other.comm()), area(other.comm()), is_glue(other.comm()) {}
    };
    
    
    template<int Dim>
    class ProjectionAlgorithm {
    public:
        using Trafo     = moonolith::AffineTransform<double, Dim-1, Dim>;
        using Shape     = moonolith::Shape<double, Dim-1, Dim>;
        using SubVector = moonolith::Vector<double, Dim-1>;
        using Vector    = moonolith::Vector<double, Dim>;
        
        ContactDataBuffers &data;
        const moonolith::IsGlue &is_glue;
        
        //algorithms
        moonolith::AffineContact<double, Dim> affine_contact;
        moonolith::WarpedContact<double, Dim> warped_contact;
        
        //buffers
        std::shared_ptr<Trafo> trafo_m, trafo_s;
        std::shared_ptr<Shape> shape_m, shape_s;
        SurfaceQuadratureConverter<Dim> converter;
        std::vector<PetscInt> dofs_petsc_s, dofs_petsc_m;
        
        //libmesh-buffers
        QMortar lm_q_master, lm_q_slave;
        std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
        libMesh::DenseMatrix<libMesh::Real> b_elmat, d_elmat;
        libMesh::DenseVector<libMesh::Real> gap_vec, normal_vec;
        libMesh::DenseMatrix<libMesh::Real> biorth_weights;
        libMesh::DenseMatrix<libMesh::Real> local_trafo;
        libMesh::DenseMatrix<libMesh::Real> inv_local_trafo;
        
        double area = 0.;
        bool use_biorth = true;
        bool use_trafo = false;
        double alpha = 1./5.;
        
        void assemble_volumes(LibMeshFunctionSpaceAdapter &space,
                              UVector &volumes
                              ) const
        {
            Chrono c;
            c.start();
            
            auto &mesh        = space.mesh();
            auto dim          = mesh.mesh_dimension();
            auto fe_type      = space.libmesh_fe_type(0);
            auto n_local_dofs = mesh.n_active_local_elem();
            
            volumes = local_zeros(n_local_dofs);
            
            {
                Write<UVector> w(volumes, utopia::LOCAL);
                
                auto fe = libMesh::FEBase::build(dim, fe_type);
                libMesh::QGauss qrule(dim, fe_type.order);
                fe->attach_quadrature_rule(&qrule);
                
                auto &JxW = fe->get_JxW();
                
                for(auto it = elements_begin(mesh); it != elements_end(mesh); ++it) {
                    fe->reinit(*it);
                    
                    auto n_qp = qrule.n_points();
                    
                    assert(n_qp > 0);
                    
                    double vol = 0.;
                    for(unsigned int qp = 0; qp < n_qp; qp++) {
                        vol += JxW[qp];
                    }
                    
                    volumes.set((*it)->id(), std::abs(vol));
                }
            }
            
            c.stop();
            std::cout << "assemble volumes (time)\n" << c << std::endl;
            double sum_v = sum(volumes);
            std::cout << "sum_v: " << sum_v << std::endl;
        }
        
        ProjectionAlgorithm(ContactDataBuffers &data, const moonolith::IsGlue &is_glue)
        : data(data), lm_q_master(Dim-1), lm_q_slave(Dim-1), is_glue(is_glue)
        {
            trafo_m = std::make_shared<Trafo>();
            trafo_s = std::make_shared<Trafo>();
            affine_contact.trafo_master = trafo_m;
            affine_contact.trafo_slave  = trafo_s;
            
            if(Dim == 2) {
                //somehow the order of surface-line elements is clockwise
                warped_contact.invert_plane_dir = true;
                
                //TODO check for the affine contact too
                // affine_contact.invert_plane_dir = true;
            }
        }
        
        void assemble_biorth_weights(const libMesh::Elem &el,
                                     const int el_order,
                                     libMesh::DenseMatrix<libMesh::Real> &weights)
        {
            std::unique_ptr<libMesh::FEBase> biorth_elem = libMesh::FEBase::build(Dim-1, libMesh::Order(el_order));
            
            const int order = order_for_l2_integral(Dim-1, el, el_order, el, el_order);
            
            libMesh::QGauss qg(Dim-1, libMesh::Order(order));
            biorth_elem->attach_quadrature_rule(&qg);
            biorth_elem->reinit(&el);
            mortar_assemble_weights(*biorth_elem, weights);
        }
        
        void init_biorth(const libMesh::Elem &el,
                         const int el_order)
        {
            if(use_biorth && biorth_weights.get_values().empty()) {
                
                use_trafo = false;
                if(el_order == 2) {
                    DualBasis::build_trafo_and_weights(el.type(), el_order, alpha, local_trafo, inv_local_trafo, biorth_weights);
                    use_trafo = true;
                } else {
                    assemble_biorth_weights(el, el_order, biorth_weights);
                }
            }
        }
        
        void init_fe(int master_order, int slave_order)
        {
            if(!master_fe) {
                master_fe = libMesh::FEBase::build(Dim-1, libMesh::Order(master_order));
                slave_fe  = libMesh::FEBase::build(Dim-1, libMesh::Order(slave_order));
                
                master_fe->get_phi();
                slave_fe->get_phi();
                slave_fe->get_JxW();
            }
        }
        
        
        void local_assemble_aux(const moonolith::Vector<double, Dim> &normal,
                                const moonolith::Storage<double> &gap)
        {
            b_elmat.zero();
            d_elmat.zero();
            gap_vec.zero();
            normal_vec.zero();
            
            if(use_biorth) {
                //P1 biorth
                mortar_assemble_weighted_biorth(
                                                *master_fe,
                                                *slave_fe,
                                                biorth_weights,
                                                b_elmat
                                                );
                
                if(use_trafo) {
                    //use trafo also for trial element
                    mortar_assemble_weighted_biorth(
                                                    *slave_fe,
                                                    local_trafo,
                                                    *slave_fe,
                                                    biorth_weights,
                                                    d_elmat
                                                    );
                    
                } else {
                    mortar_assemble_weighted_biorth(
                                                    *slave_fe,
                                                    *slave_fe,
                                                    biorth_weights,
                                                    d_elmat);
                }
                
                integrate_scalar_function_weighted_biorth(
                                                          *slave_fe,
                                                          biorth_weights,
                                                          gap,
                                                          gap_vec
                                                          );
                
                l2_project_normal_weighted_biorth(
                                                  *slave_fe,
                                                  biorth_weights,
                                                  normal,
                                                  normal_vec
                                                  );
                
                // if(!is_diag(d_elmat)) {
                //     biorth_weights.print();
                //     assert(false);
                // }
                
            } else {
                mortar_assemble(*master_fe, *slave_fe, b_elmat);
                mortar_assemble(*slave_fe, *slave_fe, d_elmat);
                integrate_scalar_function(
                                          *slave_fe,
                                          gap,
                                          gap_vec
                                          );
                
                l2_project_normal(*slave_fe, normal, normal_vec);
            }

            moonolith::print(normal, std::cout);
            
        }
        
        template<class Adapter>
        void assemble(const Adapter &master,
                      const Adapter &slave,
                      const moonolith::Quadrature<double, Dim-1> &q_master,
                      const moonolith::Quadrature<double, Dim-1> &q_slave,
                      const moonolith::Vector<double, Dim> &normal,
                      const moonolith::Storage<double> &gap
                      )
        {
            ///////////////////////////////////////////////////////
            //set-up
            
            auto &m_m = master.collection();
            auto &m_s = slave.collection();
            
            auto &e_m = master.elem();
            auto &e_s = slave.elem();
            
            auto &dofs_m = master.dofs().global;
            auto &dofs_s = slave.dofs().global;
            auto global_slave_id = slave.dofs().global_id;
            
            converter.convert_master(q_master, lm_q_master);
            converter.convert_slave(q_slave, lm_q_slave);
            
            init_fe(m_m.fe_type(0).order, m_s.fe_type(0).order);
            
            master_fe->attach_quadrature_rule(&lm_q_master);
            master_fe->reinit(&e_m);
            
            slave_fe->attach_quadrature_rule(&lm_q_slave);
            slave_fe->reinit(&e_s);
            
            init_biorth(e_s, m_s.fe_type(0).order);
            
            ///////////////////////////////////////////////////////
            //assemble local element matrices
            
            local_assemble_aux(normal, gap);
            auto isect_area = std::accumulate(
                                              b_elmat.get_values().begin(),
                                              b_elmat.get_values().end(),
                                              libMesh::Real(0.0)
                                              );
            
            ///////////////////////////////////////////////////////
            //local to global
            
            data.B.insert(dofs_s, dofs_m, b_elmat);
            data.D.insert(dofs_s, dofs_s, d_elmat);
            data.gap.insert(dofs_s, gap_vec);
            data.normal.insert_tensor_product_idx(dofs_s, Dim, normal_vec);
            data.area.insert(global_slave_id, isect_area);
            
            ///////////////////////////////////////////////////////
            //check on the area
            assert(isect_area > 0.);
            area += isect_area;
        
            auto master_tag = master.tag();
            auto slave_tag  = slave.tag();            
            
            if(is_glue(master_tag, slave_tag)) {
                data.is_glue.insert(dofs_s, 1.0);
            }

            // std::cout << isect_area << std::endl;
        }
        
        template<class Adapter>
        bool apply(Adapter &master, Adapter &slave)
        {
            auto &m_m = master.collection();
            auto &m_s = slave.collection();
            
            auto &e_m = master.elem();
            auto &e_s = slave.elem();
                
            //FIXME
            // const bool is_affine = e_m.has_affine_map() && e_s.has_affine_map();
            
            //force usage of non-affine code
            const bool is_affine = false;
            
            if(is_affine) {
                //AFFINE CONTACT
                make(e_m, affine_contact.master);
                make(e_s, affine_contact.slave);
                
                make_transform(e_m, *trafo_m);
                make_transform(e_s, *trafo_s);
                
                converter.init(
                               e_m,
                               m_m.fe_type(0).order,
                               e_s,
                               m_s.fe_type(0).order,
                               affine_contact.q_rule
                               );
                
                if(affine_contact.compute()) {
                    auto slave_area = moonolith::measure(affine_contact.slave);
                    
                    assemble(
                             master,
                             slave,
                             affine_contact.q_master,
                             affine_contact.q_slave,
                             affine_contact.normal(),
                             affine_contact.gap
                             );
                    
                    return true;
                } else {
                    return false;
                }
                
            } else {
                //WARPED CONTACT
                bool shift_in_ref_el = true;
                warped_contact.shape_master = make_shape<Dim>(e_m, m_m.libmesh_fe_type(0));
                warped_contact.shape_slave  = make_shape<Dim>(e_s, m_s.libmesh_fe_type(0));
                
                make_non_affine(e_m, warped_contact.master);
                make_non_affine(e_s, warped_contact.slave);
                
                converter.init(
                               e_m,
                               m_m.fe_type(0).order,
                               e_s,
                               m_s.fe_type(0).order,
                               warped_contact.q_rule,
                               shift_in_ref_el
                               );
                
                if(warped_contact.compute()) {
                    
                    assemble(
                             master,
                             slave,
                             warped_contact.q_master,
                             warped_contact.q_slave,
                             warped_contact.normal(),
                             warped_contact.gap
                             );
                    
                    return true;
                } else {
                    return false;
                }
            }
            
            return false;
        }
    };
    
    template<int Dim>
    class ContactAlgorithm {
    public:
        
        static bool apply(const ContactParams &params,
                          const std::shared_ptr<ElementBlackList> &black_list,
                          LibMeshFunctionSpaceAdapter &adapter,
                          ContactDataBuffers &contact_data)
        {
            using AlogrithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, LibMeshFunctionSpaceAdapter>;
            using Adapter    = typename AlogrithmT::Adapter;
            
            auto cm = std::make_shared<LibMeshCollectionManagerT>(adapter.comm(), black_list);
            
            moonolith::Communicator m_comm(cm->comm.get());
            
            moonolith::SearchSettings s;
            // s.verbosity_level = 3;
            // s.disable_redistribution = true;
            AlogrithmT algo(m_comm, cm, s);

            if(params.side_set_search_radius) {
                params.side_set_search_radius->describe(std::cout);
            
                algo.init(adapter,
                      params.contact_pair_tags,
                      params.side_set_search_radius
                    );
            } else {
                algo.init(adapter,
                      params.contact_pair_tags,
                      params.search_radius
                );
            }

            ProjectionAlgorithm<Dim> contact_algo(contact_data, *params.is_glue);
            contact_algo.use_biorth = params.use_biorthogonal_basis;
            algo.compute([&](const Adapter &master, const Adapter &slave) -> bool {
                return contact_algo.apply(master, slave);
            });
            
            m_comm.all_reduce(&contact_algo.area, 1, moonolith::MPISum());
            std::cout << "area: " << contact_algo.area << std::endl;
            
            UVector volumes;
            contact_algo.assemble_volumes(adapter,volumes);
            
            contact_data.finalize(adapter, volumes, contact_algo.alpha, contact_algo.use_biorth);
            return contact_algo.area > 0.;
        }
    };
    
    bool ContactAssembler::assemble(libMesh::MeshBase &mesh,
                                    libMesh::DofMap &dof_map,
                                    const ContactParams &params)
    {
        
        Chrono overall_time;
        overall_time.start();

        LibMeshFunctionSpaceAdapter adapter;
        
        auto spatial_dim = mesh.spatial_dimension();
        
        bool is_volume = spatial_dim == mesh.mesh_dimension();
        
        if(is_volume) {
            
            Chrono c;
            c.start();
            
            adapter.extract_surface_init_for_contact(
                                                     make_ref(mesh),
                                                     dof_map,
                                                     params.variable_number
                                                     );

            if(black_list_) {
                black_list_->init(mesh);
            }
            
            c.stop();
            
            std::cout << "adapter extract surf: " << c << std::endl;
        } else {
            //shell mesh
            assert(false); //TODO
            // adapter.init(
            //     make_ref(mesh),
            //     dof_map,
            //     params.variable_number
            // );
        }
        
        adapter.print_tags();
        ContactDataBuffers contact_data(mesh.comm().get());
        
        has_contact_ = false;
        if(spatial_dim == 2) {
            has_contact_ = ContactAlgorithm<2>::apply(params, black_list_, adapter, contact_data);
        } else if(spatial_dim == 3) {
            has_contact_ = ContactAlgorithm<3>::apply(params, black_list_, adapter, contact_data);
        }
        
        if(has_contact_) {
            contact_tensors_ = contact_data.dof_wise;
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
        std::cout << "ContactAssembler::assemble: " << overall_time << std::endl;

        return has_contact_;
    }
    
    void ContactAssembler::couple(const UVector &in, UVector &out) const
    {
        assert(contact_tensors_);
        out = transpose(contact_tensors_->complete_transformation) * in;
    }
    
    void ContactAssembler::uncouple(const UVector &in, UVector &out) const
    {
        assert(contact_tensors_);
        out = contact_tensors_->complete_transformation * in;
    }
    
    void ContactAssembler::couple(const USparseMatrix &in, USparseMatrix &out) const
    {
        assert(contact_tensors_);
        
        const auto &T = contact_tensors_->complete_transformation;
        
        out = transpose(T) * in * T;
    }
    
    const UVector &ContactAssembler::gap() const
    {
        assert(contact_tensors_);
        return contact_tensors_->gap;
    }

    UVector &ContactAssembler::gap()
    {
        assert(contact_tensors_);
        return contact_tensors_->gap;
    }

    bool ContactAssembler::init_no_contact(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map)
    {

        if(!contact_tensors_) {
            contact_tensors_ = std::make_shared<ContactTensors>();
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

    void ContactAssembler::remove_mass(const UVector &in, UVector &out) const
    {
        if(!empty(contact_tensors_->inv_mass_vector )) {
            //P1 stuff
            if(!empty(contact_tensors_->Q_inv)) {
                out = e_mul(
                    contact_tensors_->inv_mass_vector,
                    contact_tensors_->Q_inv * in
                );

                return;
            }

            out = e_mul(contact_tensors_->inv_mass_vector, in);

        } else {
            out = contact_tensors_->D_inv * in;
        }
    }

    void ContactAssembler::read(Input &in)
    {
        in.get("black-list", [this](Input &in) {
            black_list_ = std::make_shared<ElementBlackList>(true);
            black_list_->read(in);
        });
    }
}


#include "utopia_RefactoredContactTest.hpp"
#include "utopia_DualBasis.hpp"

#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_Newmark.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_ContactStabilizedNewmark.hpp"
#include "utopia_ui.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIContactParams.hpp"
#include "utopia_UIMaterial.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_InputParameters.hpp"
#include "utopia_polymorphic_QPSolver.hpp"
#include "utopia_MatrixInserter.hpp"
#include "utopia_ZeroRowsToIdentity.hpp"

#include "libmesh/mesh_refinement.h"

#include "utopia_ContactAssembler.hpp"
#include "utopia_LibMeshShape.hpp"
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

namespace utopia {

    class ContactData {
    public:
        MatrixInserter B, D;
        MatrixInserter gap, normal;
        MatrixInserter area;

        class Tensors {
        public:
            USparseMatrix B, D, Q, T;
            UVector weighted_gap, gap, normal;
            UVector area;
            UVector is_contact;

            Factorization<USparseMatrix, UVector> solver;

            void convert(
                const USparseMatrix &perm,
                const USparseMatrix &vector_perm,
                Tensors &out) const
            {
                out.B = perm * B * transpose(perm);
                out.D = perm * D * transpose(perm);

                out.weighted_gap = perm * weighted_gap;
                out.normal     = vector_perm * normal;
                out.is_contact = perm * is_contact;

                SizeType n_local_contact_nodes = 0;
                each_transform(out.is_contact, out.is_contact, [&n_local_contact_nodes](const SizeType i, const double val) -> double {
                    if(val > 0) {
                        ++n_local_contact_nodes;
                        return 1.;
                    } else {
                        return 0.0;
                    }
                });

                double sum_normal   = sum(out.normal);
                double sum_normal_e = sum(normal);

                std::cout << "n_local_contact_nodes: " << n_local_contact_nodes << std::endl;
                std::cout << "sum(normal): " << sum_normal << " == " << sum_normal_e << std::endl;

                write("P.m",  perm);
                write("Pv.m", vector_perm);
            }

            void finalize(const SizeType spatial_dim)
            {
                // D += 0. * local_identity(local_size(D));
                zero_rows_to_identity(D, 1e-13);

         
                
                gap = local_zeros(local_size(weighted_gap));
                solver.update(make_ref(D)); 
                solver.apply(weighted_gap, gap);

                // switch_dim(diag_inv_D_x, spatial_dim, diag_inv_D_xyz);
            }

            void switch_dim(
                const UVector &in,
                const SizeType spatial_dim,
                UVector &out)
            {
                out = local_zeros(local_size(in).get(0) * spatial_dim);
                Write<UVector> w(out);

                each_read(in, [&](const SizeType i, const double val) {
                    for(SizeType d = 0; d < spatial_dim; ++d) {
                        out.set(i * spatial_dim + d, val);
                    }
                });
            }
        };

        Tensors element_wise;
        Tensors dof_wise;

        inline MPI_Comm comm() const
        {
            return B.comm.get();
        }

        void finalize(
            const LibMeshFunctionSpaceAdapter &adapter,
            const UVector &volumes)
        {
            const SizeType n_local_elems = adapter.n_local_elems();
            const SizeType n_local_dofs  = adapter.n_local_dofs();
            const SizeType spatial_dim   = adapter.spatial_dim();

            area.finalize(n_local_elems);
            area.fill(element_wise.area);


            auto r = range(volumes);

            // std::vector<bool> remove(r.extent(), false);
            std::vector<bool> remove(adapter.n_local_dofs(), false);
            auto cr = adapter.permutation()->implementation().col_range();

            element_wise.is_contact = local_zeros(adapter.n_local_dofs());

            {
                Read<UVector> rv(volumes), ra(element_wise.area);
                Write<UVector> wic(element_wise.is_contact);

                for(auto i = r.begin(); i < r.end(); ++i) {
                    const auto a = element_wise.area.get(i);

                    if(a > 0.0) {
                        if(!approxeq(volumes.get(i), a, 1e-3)) {
                            remove[i - r.begin()] = true;

                            const auto &dofs = adapter.element_dof_map()[i - r.begin()].global;
                            // zeros.resize(dofs.size() * dofs.size(), 0.);
                            // element_wise.D.set_matrix(dofs, dofs, zeros);

                            for(auto d : dofs) {
                                // element_wise.weighted_gap.set(d, 0.);
                                element_wise.is_contact.set(d, 1.0);
                                remove[d - cr.begin()] = true;
                            }

                        } else {
                            std::cout << "=====================================\n";
                            std::cout << i << ") " << volumes.get(i) << " == " << element_wise.area.get(i) << std::endl;
                        }
                    }

                } 
            }


            B.finalize(n_local_dofs, n_local_dofs);
            D.finalize(n_local_dofs, n_local_dofs);
            gap.finalize(n_local_dofs);
            normal.finalize(n_local_dofs * spatial_dim);
            

            B.fill(remove, element_wise.B);
            D.fill(remove, element_wise.D);
            gap.fill(remove, element_wise.weighted_gap);
            normal.fill(element_wise.normal);


            double sum_B_x = sum(element_wise.B);
            double sum_D_x = sum(element_wise.D);
            double sum_normal_xyz = sum(element_wise.normal);

            assert(adapter.permutation());



            element_wise.convert(
                *adapter.permutation(),
                *adapter.vector_permutation(),
                dof_wise
            );

            dof_wise.finalize(spatial_dim);
        }

        ContactData(MPI_Comm comm) : B(comm), D(comm), gap(comm), normal(comm), area(comm) {}
    private:
        ContactData(const ContactData &other) : B(other.comm()), D(other.comm()), gap(other.comm()), normal(other.comm()), area(other.comm()) {}
    };


    template<int Dim>
    class SurfaceQuadratureConverter {
    public:
        using SubVector = moonolith::Vector<double, Dim-1>;
        using Vector    = moonolith::Vector<double, Dim>;

        //from libmesh to moonolith
        SubVector point_shift;
        double    point_rescale;
        double    weight_rescale;

        //from moonolith to libmesh
        double    trial_weight_rescale;
        double    test_weight_rescale;

        SubVector ref_point_shift;
        double    ref_point_rescale;

        int current_order;

        void convert_master(const moonolith::Quadrature<double, Dim-1> &in, QMortar &out)
        {
            convert(in, ref_point_shift, ref_point_rescale, trial_weight_rescale, out);
        }

        void convert_slave(const moonolith::Quadrature<double, Dim-1> &in, QMortar &out)
        {
            convert(in, ref_point_shift, ref_point_rescale, test_weight_rescale, out);
        }

        SurfaceQuadratureConverter()
        : point_rescale(1.), weight_rescale(1.), trial_weight_rescale(1.), test_weight_rescale(1.), current_order(-1)
        {}

        bool check_unity(const moonolith::Quadrature<double, Dim-1> &q)
        {
            auto sum_w = std::accumulate(q.weights.begin(), q.weights.end(), 0.);
            assert(approxeq(sum_w, 1., 1e-10));
            return approxeq(sum_w, 1., 1e-10);
        }

        void init(
            const libMesh::Elem &trial,
            const int trial_order,
            const libMesh::Elem &test,
            const int test_order,
            moonolith::Quadrature<double, Dim-1> &q,
            const bool shift_in_ref_el = false)
        {
            const int order = order_for_l2_integral(Dim-1, trial, trial_order, test, test_order);

            if(order != current_order) {
                moonolith::fill(point_shift, 0.);
                moonolith::fill(ref_point_shift, 0.);
                ref_point_rescale = 1.0;

                if(Dim == 2) {  
                    libMesh::QGauss ir(1, libMesh::Order(order));

                    if(order <= 2) {
                        ir.init(libMesh::EDGE2);
                    } else {
                        ir.init(libMesh::EDGE4);
                    }

                    point_shift.x = 1;
                    point_rescale = 0.5;
                    weight_rescale = 0.5;

                    if(shift_in_ref_el) {
                        ref_point_shift.x = -1.0;
                        ref_point_rescale = 2.0;
                    }

                    convert(ir, point_shift, point_rescale, weight_rescale, q);
                } else if(Dim == 3) {

                    libMesh::QGauss ir(2, libMesh::Order(order));

                    if(order <= 2) {
                        ir.init(libMesh::TRI3);
                    } else {
                        ir.init(libMesh::TRI6);
                    }

                    weight_rescale = 2.0;


                    if(shift_in_ref_el && is_quad(trial.type())) {
                        ref_point_shift.x = -1.0;
                        ref_point_rescale = 2.0;
                    }

                    convert(ir, point_shift, point_rescale, weight_rescale, q);

                } else {
                    assert(false);
                }

                current_order = order;
            }

            trial_weight_rescale = ref_volume(trial.type());
            test_weight_rescale  = ref_volume(test.type());

            assert(check_unity(q));
        }

    };

    template<int Dim>
    class ProjectionAlgorithm {
    public:
        using Trafo     = moonolith::AffineTransform<double, Dim-1, Dim>;
        using Shape     = moonolith::Shape<double, Dim-1, Dim>;
        using SubVector = moonolith::Vector<double, Dim-1>;
        using Vector    = moonolith::Vector<double, Dim>;

        ContactData &data;

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

        double area = 0.;
        bool use_biorth = false;
        double alpha = 1./5.;

        void assemble_volumes(
            LibMeshFunctionSpaceAdapter &space,
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

        ProjectionAlgorithm(ContactData &data)
        : data(data), lm_q_master(Dim-1), lm_q_slave(Dim-1)
        {
            trafo_m = std::make_shared<Trafo>();
            trafo_s = std::make_shared<Trafo>();
            affine_contact.trafo_master = trafo_m;
            affine_contact.trafo_slave  = trafo_s;

            if(Dim == 2) {
                //somehow the order of surface-line elements is clockwise
                warped_contact.invert_plane_dir = true;

                //TODO check for the affine contact too
            }
        }

        void assemble_biorth_weights(
                const libMesh::Elem &el,
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

        void init_biorth(
            const libMesh::Elem &el,
            const int el_order)
        {
            if(use_biorth && biorth_weights.get_values().empty()) {
                assemble_biorth_weights(el, el_order, biorth_weights);
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

        void local_assemble_aux(
            const moonolith::Vector<double, Dim> &normal,
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

                mortar_assemble_weighted_biorth(
                    *slave_fe,
                    *slave_fe, 
                    biorth_weights,
                    d_elmat);
                
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

        }

        template<class Adapter>
        void assemble(
            const Adapter &master,
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


            std::cout << moonolith::measure(warped_contact.slave) << " == " << isect_area << " == " << moonolith::measure(q_slave) << std::endl;
        }

        template<class Adapter>
        bool apply(Adapter &master, Adapter &slave)
        {
            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            auto &e_m = master.elem();
            auto &e_s = slave.elem();

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
                // bool shift_in_ref_el = false;
                // bool use_newton = false;
                // warped_contact.shape_master = std::make_shared<LibMeshShape<double, Dim>>(e_m, m_m.libmesh_fe_type(0), use_newton);
                // warped_contact.shape_slave  = std::make_shared<LibMeshShape<double, Dim>>(e_s, m_s.libmesh_fe_type(0), use_newton);

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
    bool run_contact(
        const ContactParams &params,
        LibMeshFunctionSpaceAdapter &adapter,
        ContactData &contact_data)
    {
        using AlogrithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, LibMeshFunctionSpaceAdapter>;
        using Adapter    = typename AlogrithmT::Adapter;

        auto cm = std::make_shared<LibMeshCollectionManagerT>(adapter.comm()); 

        moonolith::Communicator m_comm(cm->comm.get());

        moonolith::SearchSettings s;
        // s.verbosity_level = 3;
        // s.disable_redistribution = true;
        AlogrithmT algo(m_comm, cm, s);

        algo.init(
            adapter,
            params.contact_pair_tags,
            params.search_radius
        );

        ProjectionAlgorithm<Dim> contact_algo(contact_data);
        algo.compute([&](const Adapter &master, const Adapter &slave) -> bool {
            return contact_algo.apply(master, slave);
        });

        m_comm.all_reduce(&contact_algo.area, 1, moonolith::MPISum());
        std::cout << "area: " << contact_algo.area << std::endl;

        UVector volumes;
        contact_algo.assemble_volumes(
            adapter,
            volumes
        );

        // adapter.make_tensor_product_permutation(adapter.spatial_dim());
        contact_data.finalize(adapter, volumes);
        return contact_algo.area > 0.;
    }

    void RefactoredContactTest::run(Input &in) {

        using ProductSpaceT    = utopia::ProductFunctionSpace<LibMeshFunctionSpace>;
        using MaterialT        = utopia::UIMaterial<ProductSpaceT, USparseMatrix, UVector>;
        using ForcingFunctionT = UIForcingFunction<ProductSpaceT, UVector>;

        in.get("contact-problem", [&](Input &in) { 
            UIMesh<libMesh::DistributedMesh> mesh(this->comm());
            UIFunctionSpace<LibMeshFunctionSpace> space(make_ref(mesh));
            UIContactParams params;
            std::shared_ptr< ElasticMaterial<USparseMatrix, UVector> > model;

            in.get("mesh", mesh);
            in.get("space", space);
            in.get("contact", params);

            model = make_unique<MaterialT>(space.space());

            auto &V = space.space().subspace(0);


            LibMeshFunctionSpaceAdapter adapter;

            auto spatial_dim = V.mesh().spatial_dimension();

            bool is_volume = spatial_dim == V.mesh().mesh_dimension();

            if(is_volume) {
           
                adapter.extract_surface_init_for_contact(
                    make_ref(V.mesh()),
                    V.dof_map(),
                    params.contact_params.variable_number
                );
            } else {
                //shell mesh
                assert(false); //TODO
                // adapter.init(
                //     make_ref(V.mesh()),
                //     V.dof_map(),
                //     params.contact_params.variable_number
                // );
            }

            adapter.print_tags();

            ContactData contact_data(V.mesh().comm().get());



            bool found_contact = false;
            if(spatial_dim == 2) {
                found_contact = run_contact<2>(params.contact_params, adapter, contact_data);
            } 
            else if(spatial_dim == 3) {
                found_contact = run_contact<3>(params.contact_params, adapter, contact_data);
            }

            assert(found_contact);

            if(found_contact) {
                // write("warped.e", V, contact_data.dof_wise.gap);
                write("warped.e", V, contact_data.dof_wise.normal);
            }

            // libMesh::DenseMatrix<libMesh::Real> trafo, inv_trafo, weights;

            // DualBasis::assemble_local_trafo(libMesh::TRI6, 1./5., trafo, inv_trafo);

            // std::cout << "-----------------------\n";
            // trafo.print(std::cout);
            // std::cout << "-----------------------\n";
            // inv_trafo.print(std::cout);


            // DualBasis::assemble_biorth_weights(
            //         *V.mesh().elem(0),
            //         2,
            //         trafo,
            //         weights);

            // std::cout << "-----------------------\n";
            // weights.print(std::cout);
            // std::cout << "-----------------------\n";

        });
        
    }
}


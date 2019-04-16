#include "utopia_RefactoredContactTest.hpp"


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

#include "libmesh/mesh_refinement.h"

#include "utopia_ContactAssembler.hpp"
#include "utopia_LibMeshShape.hpp"
#include "moonolith_affine_transform.hpp"
#include "moonolith_contact.hpp"
#include "moonolith_sparse_matrix.hpp"
#include "moonolith_redistribute.hpp"

#include <vector>
#include <memory>

namespace utopia {


    //matrix proxy for utopia
    class MatrixInserter {
    public:
        MatrixInserter(MPI_Comm mpi_comm) :
          comm(mpi_comm),
          m_matrix(comm),
          redist(comm)
        {}

        template<class InsertMode>
        void finalize(InsertMode mode)
        {
            // finalize_local_structure
        }

        void fill(USparseMatrix &mat)
        {
            //
        }

        void fill(UVector &vec)
        {
            //
        }

        moonolith::Communicator comm;
        moonolith::SparseMatrix<double> m_matrix;
        moonolith::Redistribute< moonolith::SparseMatrix<double> > redist;
        std::vector<moonolith::Integer> ownership_ranges;
    };

    class ContactData {
    public:
        UVector is_contact;
        USparseMatrix B, D;

        void init(const libMesh::dof_id_type n_local_dofs)
        {
            is_contact = local_zeros(n_local_dofs);
        }

        ContactData() {}
    private:
        ContactData(const ContactData &other) {}
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

        int current_order;

        SurfaceQuadratureConverter()
        : current_order(-1)
        {}

        void init(
            const libMesh::Elem &trial,
            const int trial_order,
            const libMesh::Elem &test,
            const int test_order,
            moonolith::Quadrature<double, Dim-1> &q)
        {
            const int order = order_for_l2_integral(Dim-1, trial, trial_order, test, test_order);

            if(order != current_order) {
                moonolith::fill(point_shift, 0.);

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

                    convert(ir, point_shift, point_rescale, weight_rescale, q);
                } else if(Dim == 3) {

                    libMesh::QGauss ir(2, libMesh::Order(order));

                    if(order <= 2) {
                        ir.init(libMesh::TRI3);
                    } else {
                        ir.init(libMesh::TRI6);
                    }

                    weight_rescale = 2.0;

                    convert(ir, point_shift, point_rescale, weight_rescale, q);

                } else {
                    assert(false);
                }

                current_order = order;
            }

            trial_weight_rescale = ref_volume(trial.type());
            test_weight_rescale  = ref_volume(test.type());
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

        double area = 0.;

        ProjectionAlgorithm(ContactData &data)
        : data(data), lm_q_master(Dim-1), lm_q_slave(Dim-1)
        {
            trafo_m = std::make_shared<Trafo>();
            trafo_s = std::make_shared<Trafo>();
            affine_contact.trafo_master = trafo_m;
            affine_contact.trafo_slave  = trafo_s;


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
            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            auto &e_m = master.elem();
            auto &e_s = slave.elem();

            auto &dofs_m = master.dofs();
            auto &dofs_s = slave.dofs();

            dofs_petsc_s.clear();
            dofs_petsc_s.insert(dofs_petsc_s.end(), dofs_s.global.begin(), dofs_s.global.end());

            std::vector<double> values_s(dofs_s.global.size(), 1);
            data.is_contact.set(dofs_petsc_s, values_s);

            convert(q_master, converter.trial_weight_rescale, lm_q_master);
            convert(q_slave,  converter.test_weight_rescale, lm_q_slave);

            init_fe(m_m.fe_type(0).order, m_s.fe_type(0).order);

            master_fe->attach_quadrature_rule(&lm_q_master);
            master_fe->reinit(&e_m);

            slave_fe->attach_quadrature_rule(&lm_q_slave);
            slave_fe->reinit(&e_s);

            b_elmat.zero();
            mortar_assemble(*master_fe, *slave_fe, b_elmat);

            d_elmat.zero();
            mortar_assemble(*slave_fe, *slave_fe, d_elmat);

            auto partial_sum = std::accumulate(
                b_elmat.get_values().begin(),
                b_elmat.get_values().end(),
                libMesh::Real(0.0)
            );

            assert(partial_sum > 0.);
            area += partial_sum;
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
                bool use_newton = false;
                auto libmesh_shape = std::make_shared<LibMeshShape<double, Dim>>(e_m, m_m.libmesh_fe_type(0), use_newton);
                // libmesh_shape->verbose(true);
                warped_contact.shape_master = libmesh_shape;
               
                warped_contact.shape_slave = std::make_shared<LibMeshShape<double, Dim>>(e_s, m_s.libmesh_fe_type(0), use_newton);

                make_non_affine(e_m, warped_contact.master);
                make_non_affine(e_s, warped_contact.slave);

                converter.init(
                   e_m,
                   m_m.fe_type(0).order,
                   e_s,
                   m_s.fe_type(0).order,
                   warped_contact.q_rule
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

        Write<UVector> w(contact_data.is_contact);

        algo.compute([&](const Adapter &master, const Adapter &slave) -> bool {
            return contact_algo.apply(master, slave);
        });

        m_comm.all_reduce(&contact_algo.area, 1, moonolith::MPISum());
        std::cout << "area: " << contact_algo.area << std::endl;

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

            bool is_volume = V.mesh().spatial_dimension() == V.mesh().mesh_dimension();

            if(is_volume) {
           
                adapter.extract_surface_init(
                    make_ref(V.mesh()),
                    V.dof_map(),
                    params.contact_params.variable_number
                );
            } else {
                //shell mesh
                adapter.init(
                    make_ref(V.mesh()),
                    V.dof_map(),
                    params.contact_params.variable_number
                );
            }

            adapter.print_tags();

            ContactData contact_data;
            contact_data.init(adapter.n_local_dofs());

            bool found_contact = false;
            if(V.mesh().spatial_dimension() == 2) {
                found_contact = run_contact<2>(params.contact_params, adapter, contact_data);
            } 
            else if(V.mesh().spatial_dimension() == 3) {
                found_contact = run_contact<3>(params.contact_params, adapter, contact_data);
            }

            assert(found_contact);

            if(is_volume) {
                UVector x = (*adapter.permutation()) * contact_data.is_contact;
                write("warped.e", V, x);
            }

           
        });
        
    }
}


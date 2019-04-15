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

namespace utopia {


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

        //algorithms
        moonolith::AffineContact<double, Dim> affine_contact;
        moonolith::WarpedContact<double, Dim> warped_contact;

        //buffers
        std::shared_ptr<Trafo> trafo_m, trafo_s;
        std::shared_ptr<Shape> shape_m, shape_s;
        SurfaceQuadratureConverter<Dim> converter;

        double area = 0.;

        ProjectionAlgorithm()
        {
            trafo_m = std::make_shared<Trafo>();
            trafo_s = std::make_shared<Trafo>();
            affine_contact.trafo_master = trafo_m;
            affine_contact.trafo_slave  = trafo_s;
        }

        template<class Adapter>
        bool apply(Adapter &master, Adapter &slave)
        {
            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            auto &e_m = master.elem();
            auto &e_s = slave.elem();

            const bool is_affine = e_m.has_affine_map() && e_s.has_affine_map();

            //force usage of non-affine code
            // const bool is_affine = false;

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
                    
                    auto sum_w = 0.;
                    for(auto w : affine_contact.q_slave.weights) {
                        sum_w += w;
                    }

                    auto isect_area = sum_w * slave_area;
                    area += isect_area;
                    return true;
                } else {
                    return false;
                }

            } else {
                //WARPED CONTACT
                warped_contact.shape_master = std::make_shared<LibMeshShape<double, Dim>>(e_m, m_m.libmesh_fe_type(0));
                warped_contact.shape_slave  = std::make_shared<LibMeshShape<double, Dim>>(e_s, m_s.libmesh_fe_type(0));

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
                    auto slave_area = moonolith::measure(warped_contact.slave);
                    
                    auto sum_w = 0.;
                    for(auto w : warped_contact.q_slave.weights) {
                        sum_w += w;
                    }

                    auto isect_area = sum_w * slave_area;

                    std::cout << warped_contact.gap[0] << " " << slave_area << " " << sum_w << std::endl;
                    area += isect_area;
                    return true;
                } else {
                    return false;
                }
            }

            return false;
        }
    };

    template<int Dim>
    bool run_contact(const ContactParams &params, LibMeshFunctionSpaceAdapter &adapter)
    {
        using AlogrithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, LibMeshFunctionSpaceAdapter>;
        using Adapter    = typename AlogrithmT::Adapter;

        auto cm = std::make_shared<LibMeshCollectionManagerT>(adapter.comm()); 

        moonolith::Communicator m_comm(cm->comm.get());

        moonolith::SearchSettings s;
        // s.verbosity_level = 3;
        AlogrithmT algo(m_comm, cm, s);

        algo.init(
            adapter,
            params.contact_pair_tags,
            params.search_radius
        );

        ProjectionAlgorithm<Dim> contact_algo;

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

            if(V.mesh().spatial_dimension() == V.mesh().mesh_dimension()) {
           
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

            bool found_contact = false;
            if(V.mesh().spatial_dimension() == 2) {
                found_contact = run_contact<2>(params.contact_params, adapter);
            } 
            else if(V.mesh().spatial_dimension() == 3) {
                found_contact = run_contact<3>(params.contact_params, adapter);
            }

            assert(found_contact);
        });
        
    }
}


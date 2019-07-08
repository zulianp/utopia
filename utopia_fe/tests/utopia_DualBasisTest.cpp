#include "utopia_DualBasisTest.hpp"
#include "utopia_DualBasis.hpp"


namespace utopia {

    static void test_dual_basis(
        const libMesh::ElemType e_type,
        const int dim,
        const int order,
        const double alpha)
    {


        libMesh::DenseMatrix<libMesh::Real> trafo, inv_trafo, weights;
        DualBasis::build_trafo_and_weights(e_type, order, alpha, trafo, inv_trafo, weights);

        double sum_t = 0.0;
        
        for(auto t : trafo.get_values()) {
            sum_t += t;
        }


        auto n_nodes = trafo.m();
        // trafo.print();

        // std::cout << "------------------------\n";

        assert(approxeq(n_nodes, sum_t, 1e-10));

        double sum_w = 0.0;
        
        for(auto w : weights.get_values()) {
            sum_w += w;
        }

        assert(approxeq(n_nodes, sum_w, 1e-10));

        std::cout << "------------------------\n";
        // weights.print();

        //////////////////////////////////////////////////
        //////////////////////////////////////////////////

        libMesh::QGauss qrule(dim, libMesh::Order(2*order));
        qrule.init(e_type);

        auto fe = libMesh::FEBase::build(dim, libMesh::Order(order));
        fe->get_phi();
        fe->get_JxW();
        fe->attach_quadrature_rule(&qrule);

        auto &ref_elem = libMesh::ReferenceElem::get(e_type);
        fe->reinit(&ref_elem);

        ///////////////////////////////////////////////

        libMesh::DenseMatrix<libMesh::Real> mass_mat;
        mortar_assemble(
            *fe,
            *fe, 
            mass_mat
        );

        double sum_mm = 0.0;
        for(auto mm : mass_mat.get_values()) {
            sum_mm += mm;
        }

        // std::cout << "sum(mass_mat): " << sum_mm << std::endl;

        ///////////////////////////////////////////////
        libMesh::DenseMatrix<libMesh::Real> coupling;

        mortar_assemble_weighted_biorth(
            *fe,
            *fe, 
            weights,
            coupling
        );

        // std::cout << "------------------------\n";
        // coupling.print();

        libMesh::DenseMatrix<libMesh::Real> mass_mat_dual, inv_mass_mat_dual;
        mortar_assemble_weighted_biorth(
            *fe,
            trafo,
            *fe, 
            weights,
            mass_mat_dual
        );

        assert(is_diag(mass_mat_dual));

        double sum_mm_dual = 0.0;
        
        for(auto mm : mass_mat_dual.get_values()) {
            sum_mm_dual += mm;
        }

        // std::cout << "------------------------\n";
        // mass_mat_dual.print();

        inv_mass_mat_dual.resize(mass_mat_dual.m(), mass_mat_dual.n());
        inv_mass_mat_dual.zero();

        for(unsigned int i = 0; i < mass_mat_dual.m(); ++i) {
            for(int j = 0; j < mass_mat_dual.n(); ++j) {
                if(i == j) {
                   inv_mass_mat_dual(i, j) = 1./mass_mat_dual(i, j);
                }
            }
        }

        // std::cout << "sum(mass_mat_dual): " << sum_mm_dual << std::endl;

        assert(approxeq(sum_mm, sum_mm_dual, 1e-10));

        libMesh::DenseMatrix<libMesh::Real> trafo_t;
        trafo.get_transpose(trafo_t);

        // std::cout << "------------------------\n";
        // trafo_t.print();

        libMesh::DenseMatrix<libMesh::Real> T;

        T = coupling;
        T.left_multiply(inv_mass_mat_dual);

        // std::cout << "------------------------\n";
        // T.print();

        T.left_multiply(trafo_t);

        // std::cout << "------------------------\n";
        // T.print();

        double sum_T = 0.0;
        for(auto t : T.get_values()) {
            sum_T += t;
        }

        assert(approxeq(n_nodes, sum_T, 1e-10));
    }

    void DualBasisTest::run(Input &)
    {
        std::cout << "DualBasisTest" << std::endl;

        // libMesh::DistributedMesh mesh(comm());

        // libMesh::MeshTools::Generation::build_square(
        //     mesh,
        //     1, 1,
        //     0., 1.0,
        //     0., 1.0,
        //     libMesh::TRI6
        // );

        auto alpha = 1./5.;
        auto dim = 2;
        auto order = 2;

        // test_dual_basis(libMesh::TRI6, dim, order, alpha);
        test_dual_basis(libMesh::QUAD4, dim, 1, alpha);
        test_dual_basis(libMesh::QUAD8, dim, 2, alpha);

        dim = 1;
        test_dual_basis(libMesh::EDGE3, dim, order, alpha);

        dim = 3;
        test_dual_basis(libMesh::TET10, dim, order, alpha);

        // libMesh::DenseMatrix<libMesh::Real> trafo, inv_trafo, weights;
        // DualBasis::build_trafo_and_weights(
        //            libMesh::TET10,
        //            2,
        //            1./5,
        //            trafo,
        //            inv_trafo,
        //            weights);

        // weights.print();
        // std::cout << "------------------------\n";
        // trafo.print();



        libMesh::DenseMatrix<libMesh::Real> trafo, inv_trafo, weights;
        DualBasis::build_trafo_and_weights(
                   libMesh::HEX27,
                   2,
                   1./5,
                   trafo,
                   inv_trafo,
                   weights);

        weights.print();
        std::cout << "------------------------\n";
        trafo.print();

        std::cout << "------------------------\n";
        inv_trafo.print();
        std::cout << "------------------------\n";
        trafo.print();
    }
}


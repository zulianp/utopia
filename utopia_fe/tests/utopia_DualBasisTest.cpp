#include "utopia_DualBasisTest.hpp"
#include "utopia_DualBasis.hpp"


namespace utopia {

    void DualBasisTest::run(Input &)
    {
        std::cout << "DualBasisTest" << std::endl;

        libMesh::DistributedMesh mesh(comm());

        libMesh::MeshTools::Generation::build_square(
            mesh,
            1, 1,
            0., 1.0,
            0., 1.0,
            libMesh::TRI6
        );

        libMesh::DenseMatrix<libMesh::Real> trafo, inv_trafo, weights;
        DualBasis::build_trafo_and_weights(mesh.elem(0)->type(), 2, 1./5., trafo, inv_trafo, weights);

        double sum_t = 0.0;
        
        for(auto t : trafo.get_values()) {
            sum_t += t;
        }

        // trafo.print();

        // std::cout << "------------------------\n";

        assert(approxeq(6.0, sum_t, 1e-10));

        double sum_w = 0.0;
        
        for(auto w : weights.get_values()) {
            sum_w += w;
        }

        assert(approxeq(6.0, sum_w, 1e-10));

        // weights.print();

        //////////////////////////////////////////////////
        //////////////////////////////////////////////////


        libMesh::QGauss qrule(2, libMesh::Order(4));
        qrule.init(libMesh::TRI6);

        auto fe = libMesh::FEBase::build(2, libMesh::Order(2));
        fe->get_phi();
        fe->get_JxW();
        fe->attach_quadrature_rule(&qrule);

        auto &ref_elem = libMesh::ReferenceElem::get(mesh.elem(0)->type());
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

        libMesh::DenseMatrix<libMesh::Real> mass_mat_dual;
        mortar_assemble_weighted_biorth(
            *fe,
            trafo,
            *fe, 
            weights,
            mass_mat_dual
        );

        double sum_mm_dual = 0.0;
        
        for(auto mm : mass_mat_dual.get_values()) {
            sum_mm_dual += mm;
        }

        // std::cout << "sum(mass_mat_dual): " << sum_mm_dual << std::endl;

        assert(approxeq(sum_mm, sum_mm_dual, 1e-10));
    }
}
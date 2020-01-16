#ifndef UTOPIA_SURFACE_QUADRATURE_CONVERTER_HPP
#define UTOPIA_SURFACE_QUADRATURE_CONVERTER_HPP

#include "moonolith_vector.hpp"
#include "moonolith_map_quadrature.hpp"
#include "utopia_LibMeshToMoonolithConvertions.hpp"

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
}


#endif //UTOPIA_SURFACE_QUADRATURE_CONVERTER_HPP

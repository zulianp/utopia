

#include "utopia_Main.hpp"

#include "utopia_intrepid2.hpp"
#include "utopia_stk.hpp"

#include "utopia_kokkos_SamplePoints.hpp"

namespace utopia {

    template <class FunctionSpace>
    class SamplePointsApp {
    public:
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Field_t = utopia::Field<FunctionSpace>;
        using Intrepid2FE_t = utopia::intrepid2::FE<Scalar_t>;
        using DynRankView_t = typename Intrepid2FE_t::DynRankView;
        using Intrepid2Field_t = utopia::kokkos::Field<Intrepid2FE_t>;
        using SamplePoints_t = utopia::kokkos::SamplePoints<Intrepid2FE_t>;
        using HPolytope_t = utopia::kokkos::HPolytope<DynRankView_t>;

        void run(Input &in) {
            FunctionSpace space;
            space.read(in);

            if (space.empty()) {
                space.comm().root_print("No input space and mesh provided, using unit cube!\n");

                space.mesh().unit_cube(3, 3, 3);
                space.initialize();
            }

            auto fe = std::make_shared<Intrepid2FE_t>();

            DynRankView_t normals, distances;
            HPolytope_t h_poly(normals, distances);
            SamplePoints_t sample(fe, h_poly);
            sample.read(in);

            DynRankView_t values;
            Intrepid2Field_t field(fe);
            sample.sample(field, values);
        }
    };

}  // namespace utopia

void stk_sample_points(utopia::Input &in) {
    utopia::SamplePointsApp<utopia::stk::FunctionSpace> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(stk_sample_points);

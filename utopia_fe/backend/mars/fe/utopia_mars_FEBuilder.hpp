#ifndef UTOPIA_MARS_FE_BUILDER_HPP
#define UTOPIA_MARS_FE_BUILDER_HPP

#include "mars_quad4.hpp"

#include "utopia_kokkos_UniformFE.hpp"

namespace utopia {
    namespace mars {
        template <class Handler, class FE>
        class FEBuilder {
        public:
        };

        template <class Handler, typename... Args>
        class FEBuilder<Handler, utopia::kokkos::UniformFE<Traits<mars::FunctionSpace>::Scalar, Args...>> {
        public:
            using Matrix = Traits<mars::FunctionSpace>::Matrix;
            using Vector = Traits<mars::FunctionSpace>::Vector;
            using Scalar = Traits<mars::FunctionSpace>::Scalar;

            using DofHandler = typename Handler::DofHandler;
            using FEDofMap = typename Handler::FEDofMap;
            using SPattern = typename Handler::SPattern;

            static const ::mars::Integer Type = SPattern::DofHandler::ElemType;

            static constexpr int Dim = Handler::Dim;

            using FE = utopia::kokkos::UniformFE<Scalar>;
            using UniformFE = typename utopia::mars::FETypeSelect<Scalar, Dim>::Type;

            void build(Handler &handler, FE &fe) {
                using Coordinates = typename UniformFE::Coordinates;

                Coordinates coords("coords");

                // auto sp = handler.get_sparsity_pattern();
                auto fe_dof_map = handler.get_fe_dof_map();
                auto dof_handler = handler.get_dof_handler();
                int block_size = dof_handler.get_block();

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    if (elem_index == 0) {
                        Scalar p[Dim];

                        for (int i = 0; i < coords.extent(0); ++i) {
                            auto local_dof = fe_dof_map.get_elem_local_dof(elem_index, i * block_size);
                            dof_handler.template get_dof_coordinates_from_local<Type>(local_dof, p);

                            // printf("%d ", int(local_dof));

                            for (int d = 0; d < coords.extent(1); ++d) {
                                coords(i, d) = p[d];

                                // printf("%g ", (double) p[d]);
                            }

                            // printf("\n");
                        }
                    }
                });

                auto deprecated_fe = utopia::make_unique<UniformFE>();
                deprecated_fe->init(coords);

                typename FE::MeasureView measure("measure", deprecated_fe->measure.extent(0));
                ::Kokkos::deep_copy(measure, deprecated_fe->measure);
                typename FE::FunctionView fun("fun", deprecated_fe->fun.extent(0), deprecated_fe->fun.extent(1));
                ::Kokkos::deep_copy(fun, deprecated_fe->fun);
                typename FE::GradientView grad("grad",
                                               deprecated_fe->grad.extent(0),
                                               deprecated_fe->grad.extent(1),
                                               deprecated_fe->grad.extent(2));
                ::Kokkos::deep_copy(grad, deprecated_fe->grad);

                // TODO
                typename FE::JacobianView jacobian("jacobian", 0, 0, 0);
                // ::Kokkos::deep_copy(jacobian, deprecated_fe->jacobian);
                typename FE::JacobianInverseView jacobian_inverse("jacobian_inverse", 0, 0, 0);
                // ::Kokkos::deep_copy(jacobian_inverse, deprecated_fe->jacobian_inverse);

                // auto n_elems = dof_handler.get_elem_size();
                auto n_elems = fe_dof_map.get_fe_dof_map_size();
                fe.init(n_elems, measure, fun, grad, jacobian, jacobian_inverse);
            }
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_FE_BUILDER_HPP
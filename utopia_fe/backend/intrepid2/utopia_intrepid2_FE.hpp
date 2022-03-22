#ifndef UTOPIA_INTREPID2_FE_HPP
#define UTOPIA_INTREPID2_FE_HPP

#include "utopia_fe_base.hpp"

#include "utopia_kokkos_FE.hpp"

#include "utopia_intrepid2_Base.hpp"
#include "utopia_intrepid2_ShellTools.hpp"

#include <Kokkos_DynRankView.hpp>
#include <Shards_CellTopology.hpp>

#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_Cubature.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>
#include <Intrepid2_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid2_HGRAD_QUAD_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TET_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_C1_FEM.hpp>

namespace utopia {

    namespace intrepid2 {

        template <typename Scalar, class ExecutionSpace_ = utopia::intrepid2::ExecutionSpace>
        class Intrepid2FETraits {
        public:
            using ExecutionSpace = ExecutionSpace_;
            using DynRankView = utopia::intrepid2::ViewDevice<Scalar>;
            using SizeType = int;
            using MeasureView = DynRankView;
            using FunctionView = DynRankView;
            using GradientView = DynRankView;
            using JacobianView = DynRankView;
            using JacobianInverseView = DynRankView;
            using IntView = utopia::intrepid2::IntViewDevice;
        };

        template <typename Scalar_>
        class FE : public utopia::kokkos::FE<Scalar_, Intrepid2FETraits<Scalar_>> {
        public:
            using Scalar = Scalar_;
            using Super = utopia::kokkos::FE<Scalar, Intrepid2FETraits<Scalar>>;
            using HostExecutionSpace = utopia::intrepid2::HostExecutionSpace;
            using ExecutionSpace = utopia::intrepid2::ExecutionSpace;
            using DynRankView = utopia::intrepid2::ViewDevice<Scalar>;
            using IntView = utopia::intrepid2::IntViewDevice;

            using Cubature = ::Intrepid2::Cubature<ExecutionSpace, Scalar, Scalar>;
            using CubaturePtr = ::Teuchos::RCP<Cubature>;
            using CellTopology = ::shards::CellTopology;
            using Tet = ::Intrepid2::Basis_HGRAD_TET_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using Tri = ::Intrepid2::Basis_HGRAD_TRI_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using Quad = ::Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using Line = ::Intrepid2::Basis_HGRAD_LINE_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using Hex = ::Intrepid2::Basis_HGRAD_HEX_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using CellTools = ::Intrepid2::CellTools<ExecutionSpace>;
            using FunctionSpaceTools = ::Intrepid2::FunctionSpaceTools<ExecutionSpace>;
            using SizeType = ::Intrepid2::ordinal_type;

            using CellTestTrialRange = Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>;
            using CellTestRange = Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>;
            using CellRange = Kokkos::RangePolicy<ExecutionSpace>;
            using CellQPRange = Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>;

            virtual ~FE() = default;

            void init(const CellTopology &type, const DynRankView &cell_nodes, const int degree = 1) {
                this->type = type;
                this->cell_nodes = cell_nodes;

                // ::Intrepid2::DefaultCubatureFactory<Scalar> cub_factory;
                ::Intrepid2::DefaultCubatureFactory cub_factory;
                cubature = cub_factory.create<ExecutionSpace, Scalar, Scalar>(type, degree);

                switch (type.getKey()) {
                    case shards::Line<>::key: {
                        Line line;
                        init_aux(line);
                        break;
                    }

                    case shards::Triangle<>::key: {
                        Tri tri;
                        init_aux(tri);
                        break;
                    }

                    case shards::Quadrilateral<>::key: {
                        Quad quad;
                        init_aux(quad);
                        break;
                    }

                    case shards::Tetrahedron<>::key: {
                        Tet tet;
                        init_aux(tet);
                        break;
                    }

                    case shards::Hexahedron<>::key: {
                        Hex hex;
                        init_aux(hex);
                        break;
                    }

                    default: {
                        utopia::err() << "Unsupported type " << type << '\n';
                        assert(false);
                        Utopia::Abort("FE init failed!");
                        break;
                    }
                }
            }

            inline SizeType n_cells() const override { return cell_nodes.extent(0); }
            // inline SizeType n_shape_functions() const override { return grad.extent(1); }
            inline SizeType n_quad_points() const override { return cubature->getNumPoints(); }
            inline SizeType spatial_dimension() const override { return cell_nodes.extent(2); }
            inline SizeType manifold_dimension() const override { return type.getDimension(); }
            inline DynRankView points() const { return cell_nodes; }

            CellTopology type;
            DynRankView cell_nodes;
            CubaturePtr cubature;
            DynRankView q_points;
            DynRankView q_weights;
            DynRankView jacobian_det;

            // NVCC_PRIVATE :
            template <class BasisType>
            void init_aux(const BasisType &basis) {
                if (this->is_shell()) {
                    init_aux_shell(basis);
                    return;
                }

                DynRankView fun;
                DynRankView grad;
                DynRankView jacobian;
                DynRankView jacobian_inv;
                DynRankView measure;

                auto spatial_dimension = this->spatial_dimension();
                auto manifold_dimension = this->manifold_dimension();

                auto n_quad_points = this->n_quad_points();
                auto n_cells = this->n_cells();
                auto n_fun = basis.getCardinality();

                assert(n_cells > 0);
                assert(n_quad_points > 0);
                assert(spatial_dimension > 0);
                assert(manifold_dimension > 0);

                q_weights = DynRankView("q_weights", n_quad_points);
                q_points = DynRankView("q_points", n_quad_points, manifold_dimension);
                jacobian = DynRankView("jacobian", n_cells, n_quad_points, spatial_dimension, spatial_dimension);
                jacobian_inv =
                    DynRankView("jacobian_inv", n_cells, n_quad_points, spatial_dimension, spatial_dimension);
                measure = DynRankView("measure", n_cells, n_quad_points);
                jacobian_det = DynRankView("jacobian_det", n_cells, n_quad_points);

                cubature->getCubature(q_points, q_weights);
                CellTools::setJacobian(jacobian, q_points, cell_nodes, type);
                CellTools::setJacobianInv(jacobian_inv, jacobian);
                CellTools::setJacobianDet(jacobian_det, jacobian);

                FunctionSpaceTools::computeCellMeasure<Scalar>(measure, jacobian_det, q_weights);

                fun = DynRankView("fun", n_fun, n_quad_points);
                basis.getValues(fun, q_points, ::Intrepid2::OPERATOR_VALUE);

                DynRankView ref_grad("ref_grad", n_fun, n_quad_points, manifold_dimension);
                grad = DynRankView("grad", n_cells, n_fun, n_quad_points, spatial_dimension);

                basis.getValues(ref_grad, q_points, ::Intrepid2::OPERATOR_GRAD);
                FunctionSpaceTools::HGRADtransformGRAD<Scalar>(grad, jacobian_inv, ref_grad);

                Super::init(measure, fun, grad, jacobian, jacobian_inv);
            }

            template <class BasisType>
            void init_aux_shell(const BasisType &basis) {
                DynRankView fun;
                DynRankView grad;
                DynRankView jacobian;
                DynRankView jacobian_inv;
                DynRankView measure;

                auto spatial_dimension = this->spatial_dimension();
                auto manifold_dimension = this->manifold_dimension();

                auto n_quad_points = this->n_quad_points();
                auto n_cells = this->n_cells();
                auto n_fun = basis.getCardinality();

                assert(n_cells > 0);
                assert(n_quad_points > 0);
                assert(spatial_dimension > 0);
                assert(manifold_dimension > 0);

                DynRankView q_weights("q_weights", n_quad_points);
                q_points = DynRankView("q_points", n_quad_points, manifold_dimension);
                cubature->getCubature(q_points, q_weights);

                fun = DynRankView("fun", n_fun, n_quad_points);
                basis.getValues(fun, q_points, ::Intrepid2::OPERATOR_VALUE);

                DynRankView ref_grad("ref_grad", n_fun, n_quad_points, manifold_dimension);
                grad = DynRankView("grad", n_cells, n_fun, n_quad_points, spatial_dimension);

                basis.getValues(ref_grad, q_points, ::Intrepid2::OPERATOR_GRAD);

                ShellTools<Scalar>::allocate_jacobian(
                    n_cells, manifold_dimension, spatial_dimension, n_quad_points, jacobian);
                ShellTools<Scalar>::allocate_jacobian_inverse(
                    n_cells, manifold_dimension, spatial_dimension, n_quad_points, jacobian_inv);
                ShellTools<Scalar>::allocate_measure(n_cells, n_quad_points, measure);
                ShellTools<Scalar>::cell_geometry(cell_nodes, q_weights, ref_grad, jacobian, jacobian_inv, measure);
                ShellTools<Scalar>::transform_gradient_to_physical_space(jacobian_inv, ref_grad, grad);

                Super::init(measure, fun, grad, jacobian, jacobian_inv);
                // print_function();
                // print_gradient();
                // print_measure();
                // print_jacobian();
                // print_jacobian_inverse();
            }
        };

    }  // namespace intrepid2

}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FE_HPP

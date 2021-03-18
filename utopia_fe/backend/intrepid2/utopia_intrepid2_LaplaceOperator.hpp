#ifndef UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP
#define UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP

#include "utopia_fe_base.hpp"

#include <Kokkos_DynRankView.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_Cubature.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>
#include <Intrepid2_HGRAD_TET_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_C1_FEM.hpp>

namespace utopia {

    template <class DiffusionCoefficient>
    class LaplaceOperator {
    public:
        DiffusionCoefficient coeff;
    };

    namespace intrepid2 {

        template <class Operator, typename Scalar = UScalar>
        class Assemble {};

        template <typename Scalar>
        class HGradFE {
        public:
            using ExecutionSpace = ::Kokkos::Serial;
            // using ExecutionSpace = ::Kokkos::DefaultExecutionSpace;
            using DynRankView = ::Kokkos::DynRankView<Scalar>;
            using Cubature = ::Intrepid2::Cubature<ExecutionSpace, Scalar, Scalar>;
            using CubaturePtr = ::Teuchos::RCP<Cubature>;
            using CellTopology = ::shards::CellTopology;
            using Tet = ::Intrepid2::Basis_HGRAD_TET_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using Tri = ::Intrepid2::Basis_HGRAD_TRI_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using CellTools = ::Intrepid2::CellTools<ExecutionSpace>;
            using FunctionSpaceTools = ::Intrepid2::FunctionSpaceTools<ExecutionSpace>;

            virtual ~HGradFE() = default;

            void init(const CellTopology &type, const DynRankView &cell_nodes, const int degree = 1) {
                this->type = type;
                this->cell_nodes = cell_nodes;

                // ::Intrepid2::DefaultCubatureFactory<Scalar> cub_factory;
                ::Intrepid2::DefaultCubatureFactory cub_factory;
                cubature = cub_factory.create<ExecutionSpace, Scalar, Scalar>(type, degree);

                switch (type.getKey()) {
                    case shards::Line<>::key: {
                        assert(false);
                        break;
                    }

                    case shards::Triangle<>::key: {
                        Tri tri;
                        init_aux(tri);
                        break;
                    }

                    case shards::Tetrahedron<>::key: {
                        Tet tet;
                        init_aux(tet);
                        break;
                    }

                    default: {
                        assert(false);
                        break;
                    }
                }
            }

            CellTopology type;
            DynRankView cell_nodes;
            CubaturePtr cubature;

            DynRankView grad;
            DynRankView jacobian;
            DynRankView jacobian_inv;
            DynRankView jacobian_det;
            DynRankView q_points;
            DynRankView measure;

            // NVCC_PRIVATE :
            template <class BasisType>
            void init_aux(const BasisType &basis) {
                auto spatial_dimension = type.getDimension();
                auto num_nodes = type.getNodeCount();
                auto num_qp = cubature->getNumPoints();
                auto num_cells = cell_nodes.extent(0);

                assert(num_cells > 0);
                assert(num_qp > 0);
                assert(num_nodes > 0);
                assert(spatial_dimension > 0);

                DynRankView q_weights("q_weights", num_qp);
                q_points = DynRankView("q_points", num_qp, spatial_dimension);
                jacobian = DynRankView("jacobian", num_cells, num_qp, spatial_dimension, spatial_dimension);
                jacobian_inv = DynRankView("jacobian_inv", num_cells, num_qp, spatial_dimension, spatial_dimension);
                measure = DynRankView("measure", num_cells, num_qp);
                jacobian_det = DynRankView("jacobian_det", num_cells, num_qp);

                cubature->getCubature(q_points, q_weights);
                CellTools::setJacobian(jacobian, q_points, cell_nodes, type);
                CellTools::setJacobianInv(jacobian_inv, jacobian);
                CellTools::setJacobianDet(jacobian_det, jacobian);

                FunctionSpaceTools::computeCellMeasure<Scalar>(measure, jacobian_det, q_weights);
            }
        };

        template <class DiffusionCoefficient, typename Scalar>
        class Assemble<LaplaceOperator<DiffusionCoefficient>, Scalar> {
        public:
        };

    }  // namespace intrepid2
}  // namespace utopia

#endif
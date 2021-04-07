#ifndef UTOPIA_INTREPID2_FE_HPP
#define UTOPIA_INTREPID2_FE_HPP

#include "utopia_fe_base.hpp"

#include <Kokkos_DynRankView.hpp>

// #include <KokkosBlas1_scal.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_Cubature.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>
#include <Intrepid2_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TET_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_C1_FEM.hpp>

namespace utopia {

    namespace intrepid2 {

        template <typename Scalar>
        class FE {
        public:
            using ExecutionSpace = ::Kokkos::Serial;
            // using ExecutionSpace = ::Kokkos::DefaultExecutionSpace;
            using DynRankView = ::Kokkos::DynRankView<Scalar>;
            using Cubature = ::Intrepid2::Cubature<ExecutionSpace, Scalar, Scalar>;
            using CubaturePtr = ::Teuchos::RCP<Cubature>;
            using CellTopology = ::shards::CellTopology;
            using Tet = ::Intrepid2::Basis_HGRAD_TET_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using Tri = ::Intrepid2::Basis_HGRAD_TRI_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using Hex = ::Intrepid2::Basis_HGRAD_HEX_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using CellTools = ::Intrepid2::CellTools<ExecutionSpace>;
            using FunctionSpaceTools = ::Intrepid2::FunctionSpaceTools<ExecutionSpace>;
            using SizeType = ::Intrepid2::ordinal_type;

            virtual ~FE() = default;

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

                    case shards::Hexahedron<>::key: {
                        Hex hex;
                        init_aux(hex);
                        break;
                    }

                    default: {
                        assert(false);
                        break;
                    }
                }
            }

            inline SizeType num_cells() const { return cell_nodes.extent(0); }
            inline SizeType num_fields() const { return grad.extent(1); }
            inline SizeType num_qp() const { return cubature->getNumPoints(); }
            inline SizeType spatial_dimension() const { return cell_nodes.extent(2); }
            inline SizeType manifold_dimension() const { return type.getDimension(); }

            bool is_trace() const { return manifold_dimension() < spatial_dimension(); }

            CellTopology type;
            // n_local_elements x n_nodes_x_elem x spatial_dim;
            DynRankView cell_nodes;
            CubaturePtr cubature;

            DynRankView fun;
            DynRankView grad;
            DynRankView jacobian;
            DynRankView jacobian_inv;
            DynRankView jacobian_det;
            DynRankView q_points;
            DynRankView measure;

            // NVCC_PRIVATE :
            template <class BasisType>
            void init_aux(const BasisType &basis) {
                auto spatial_dimension = this->spatial_dimension();
                auto manifold_dimension = this->manifold_dimension();

                auto num_qp = this->num_qp();
                auto num_cells = this->num_cells();
                auto n_fun = basis.getCardinality();

                assert(num_cells > 0);
                assert(num_qp > 0);
                assert(spatial_dimension > 0);
                assert(manifold_dimension > 0);

                DynRankView q_weights("q_weights", num_qp);
                q_points = DynRankView("q_points", num_qp, manifold_dimension);
                jacobian = DynRankView("jacobian", num_cells, num_qp, spatial_dimension, manifold_dimension);
                jacobian_inv = DynRankView("jacobian_inv", num_cells, num_qp, manifold_dimension, spatial_dimension);
                measure = DynRankView("measure", num_cells, num_qp);
                jacobian_det = DynRankView("jacobian_det", num_cells, num_qp);

                cubature->getCubature(q_points, q_weights);
                CellTools::setJacobian(jacobian, q_points, cell_nodes, type);
                CellTools::setJacobianInv(jacobian_inv, jacobian);
                CellTools::setJacobianDet(jacobian_det, jacobian);

                FunctionSpaceTools::computeCellMeasure<Scalar>(measure, jacobian_det, q_weights);

                fun = DynRankView("fun", n_fun, num_qp);
                basis.getValues(fun, q_points, ::Intrepid2::OPERATOR_VALUE);

                DynRankView ref_grad("ref_grad", n_fun, num_qp, manifold_dimension);
                grad = DynRankView("grad", num_cells, n_fun, num_qp, spatial_dimension);

                basis.getValues(ref_grad, q_points, ::Intrepid2::OPERATOR_GRAD);
                FunctionSpaceTools::HGRADtransformGRAD<Scalar>(grad, jacobian_inv, ref_grad);
            }
        };
    }  // namespace intrepid2

}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FE_HPP
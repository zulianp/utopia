#ifndef UTOPIA_INTREPID2_FE_HPP
#define UTOPIA_INTREPID2_FE_HPP

#include "utopia_fe_base.hpp"

#include "utopia_intrepid2_ShellTools.hpp"

#include <Kokkos_DynRankView.hpp>
#include <Shards_CellTopology.hpp>

#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_Cubature.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>
#include <Intrepid2_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TET_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_C1_FEM.hpp>

namespace utopia {

    template <typename Scalar_>
    class Traits<::Kokkos::DynRankView<Scalar_>> {
    public:
        using Scalar = Scalar_;
    };

    namespace intrepid2 {

        template <typename Scalar>
        class FE {
        public:
            // using ExecutionSpace = ::Kokkos::Serial;
            using HostExecutionSpace = ::Kokkos::DefaultHostExecutionSpace;
            using ExecutionSpace = ::Kokkos::DefaultExecutionSpace;
            using DynRankView = ::Kokkos::DynRankView<Scalar>;
            using IntView = ::Kokkos::DynRankView<int>;

            using HostDynRankView = ::Kokkos::DynRankView<Scalar, HostExecutionSpace>;
            using HostIntView = ::Kokkos::DynRankView<int, HostExecutionSpace>;

            using Cubature = ::Intrepid2::Cubature<ExecutionSpace, Scalar, Scalar>;
            using CubaturePtr = ::Teuchos::RCP<Cubature>;
            using CellTopology = ::shards::CellTopology;
            using Tet = ::Intrepid2::Basis_HGRAD_TET_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using Tri = ::Intrepid2::Basis_HGRAD_TRI_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using Line = ::Intrepid2::Basis_HGRAD_LINE_C1_FEM<ExecutionSpace, Scalar, Scalar>;
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
                        Line line;
                        init_aux(line);
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

            inline bool is_shell() const { return manifold_dimension() < spatial_dimension(); }

            void print_measure() {
                auto num_qp = this->num_qp();
                auto num_cells = this->num_cells();

                // Avoid capturing (this)
                auto measure = this->measure;

                Kokkos::parallel_for(
                    "FE::print_measure", num_cells, KOKKOS_LAMBDA(const int &cell) {
                        printf("cell: %d\n", cell);

                        for (int qp = 0; qp < num_qp; ++qp) {
                            printf("%g ", measure(cell, qp));
                        }

                        printf("\n");
                    });
            }

            void print_jacobian() {
                int num_cells = this->num_cells();
                int num_qp = this->num_qp();

                int spatial_dimension = this->spatial_dimension();
                int manifold_dimension = this->manifold_dimension();

                // Avoid capturing (this)
                auto jacobian = this->jacobian;

                Kokkos::parallel_for(
                    "FE::print_jacobian", num_cells, KOKKOS_LAMBDA(const int &cell) {
                        printf("cell: %d\n", cell);

                        for (int qp = 0; qp < num_qp; ++qp) {
                            for (int r = 0; r < spatial_dimension; ++r) {
                                for (int c = 0; c < manifold_dimension; ++c) {
                                    printf("%g ", jacobian(cell, qp, r, c));
                                }
                                printf("\n");
                            }

                            printf("\n");
                        }

                        printf("\n");
                    });
            }

            void print_jacobian_inverse() {
                int num_cells = this->num_cells();
                int num_qp = this->num_qp();

                int spatial_dimension = this->spatial_dimension();
                int manifold_dimension = this->manifold_dimension();

                // Avoid capturing (this)
                auto jacobian_inv = this->jacobian_inv;

                Kokkos::parallel_for(
                    "FE::print_jacobian_inverse", num_cells, KOKKOS_LAMBDA(const int &cell) {
                        printf("cell: %d\n", cell);

                        for (int qp = 0; qp < num_qp; ++qp) {
                            for (int r = 0; r < manifold_dimension; ++r) {
                                for (int c = 0; c < spatial_dimension; ++c) {
                                    printf("%g ", jacobian_inv(cell, qp, r, c));
                                }
                                printf("\n");
                            }

                            printf("\n");
                        }

                        printf("\n");
                    });
            }

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

            // Optional
            IntView element_tags;

            inline bool has_element_tags() const { return element_tags.size() > 0; }

            // NVCC_PRIVATE :
            template <class BasisType>
            void init_aux(const BasisType &basis) {
                if (is_shell()) {
                    init_aux_shell(basis);
                    return;
                }

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
                jacobian = DynRankView("jacobian", num_cells, num_qp, spatial_dimension, spatial_dimension);
                jacobian_inv = DynRankView("jacobian_inv", num_cells, num_qp, spatial_dimension, spatial_dimension);
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

            template <class BasisType>
            void init_aux_shell(const BasisType &basis) {
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
                cubature->getCubature(q_points, q_weights);

                fun = DynRankView("fun", n_fun, num_qp);
                basis.getValues(fun, q_points, ::Intrepid2::OPERATOR_VALUE);

                DynRankView ref_grad("ref_grad", n_fun, num_qp, manifold_dimension);
                grad = DynRankView("grad", num_cells, n_fun, num_qp, spatial_dimension);

                basis.getValues(ref_grad, q_points, ::Intrepid2::OPERATOR_GRAD);

                ShellTools<Scalar>::allocate_jacobian(num_cells, spatial_dimension, num_qp, jacobian);
                ShellTools<Scalar>::allocate_jacobian_inverse(num_cells, spatial_dimension, num_qp, jacobian_inv);
                ShellTools<Scalar>::allocate_measure(num_cells, num_qp, measure);
                ShellTools<Scalar>::cell_geometry(cell_nodes, jacobian, jacobian_inv, measure);
                ShellTools<Scalar>::transform_gradient_to_physical_space(jacobian_inv, ref_grad, grad);
            }
        };
    }  // namespace intrepid2

}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FE_HPP

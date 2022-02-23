#ifndef UTOPIA_MARS_FE_HPP
#define UTOPIA_MARS_FE_HPP

#include "mars.hpp"
namespace utopia {
    namespace mars {

        template <typename Real, class SparsityPattern, class FEM>
        class FE {
        public:
            using Integer = ::mars::Integer;
            using ViewVectorType = ::mars::ViewVectorType<Real>;
            using ViewMatrixType = ::mars::ViewMatrixType<Real>;

            // FIXME Hardcoded
            using Elem = ::mars::FEQuad4<Real>;
            // static const ::mars::Integer Type = SparsityPattern::DofHandler::ElemType;
            static const ::mars::Integer Type = 4;
            static const int Dim = 2;
            // Hardcoded

            using Quadrature = typename Elem::Quadrature;

            // use as more readable tuple index to identify the data
            FE(SparsityPattern &d, FEM &f) : sparsity_pattern_(d), fe_(f) {}

            // OK as long as we are doing affine meshes in 2D

            static void compute_invJ_and_detJ(const SparsityPattern &data,
                                              const FEM &fe,
                                              ViewVectorType detJ,
                                              ViewMatrixType invJ) {
                auto dof_handler = data.get_dof_handler();

                fe.iterate(MARS_LAMBDA(const Integer elem_index) {
                    Real J[Dim * Dim];
                    Real J_inv[Dim * Dim];
                    Real point_ref[Dim];
                    Real p2[Dim];
                    Real p3[Dim];

                    Integer local_dof_0 = fe.get_elem_local_dof(elem_index, 0);
                    dof_handler.template get_dof_coordinates_from_local<Type>(local_dof_0, point_ref);

                    Integer local_dof_1 = fe.get_elem_local_dof(elem_index, 1);
                    dof_handler.template get_dof_coordinates_from_local<Type>(local_dof_1, p2);

                    // col 0, p1
                    J[0] = p2[0] - point_ref[0];
                    J[2] = p2[1] - point_ref[1];

                    // we skip p2 and get p3
                    Integer local_dof_3 = fe.get_elem_local_dof(elem_index, 3);
                    dof_handler.template get_dof_coordinates_from_local<Type>(local_dof_3, p3);

                    // col 1, p3
                    J[1] = p3[0] - point_ref[0];
                    J[3] = p3[1] - point_ref[1];

                    // compute the determinant and the inverse of the Jacobian
                    Real det_J;
                    ::mars::Invert<Dim>::apply(J, J_inv, det_J);

                    assert(det_J > 0);

                    detJ(elem_index) = Kokkos::ArithTraits<Real>::abs(det_J);

                    for (int i = 0; i < 4; i++) {
                        invJ(elem_index, i) = J_inv[i];
                    }
                });

                // Real measure = KokkosBlas::nrm1(det_J_);
                // std::cout << "measure: " << measure << std::endl;
            }

            // template <Integer INPUT>
            // MARS_INLINE_FUNCTION void gather_elem_data(const SparsityPattern &sp,
            //                                            const FEM &fe,
            //                                            const Integer elem_index,
            //                                            DMDataType<INPUT> *sol) {
            //     for (int i = 0; i < FEM::elem_nodes; i++) {
            //         // forach dof get the local number
            //         const Integer local_dof = fe.get_elem_local_dof(elem_index, i);
            //         // use the local number to read the corresponding user data
            //         sol[i] = sp.template get_dof_data<INPUT>(local_dof);
            //     }
            // }

            // if non-linear than the quad rule should be computed per quad point and not
            // anymore for each element so the coalescing of the Jinv will not matter.
            // In that case maybe a better way to go is parallel through the quad points.

            // void intrepid2_laplace_operator() {
            //     integrate_laplace_operator(sparsity_pattern_, fe_, quad, det_J_, inv_J_);
            // }

            // void integrate_laplace_operator(const SparsityPattern &sp,
            //                                 const FEM &fe,
            //                                 const Quadrature &quad,
            //                                 ViewVectorType det_J,
            //                                 ViewMatrixType J_inv) {
            //     constexpr int n_qp = Quadrature::n_points();
            //     constexpr int dim = Quadrature::dim();

            //     ::mars::ViewVectorTextureC<Real, n_qp> q_weights = quad.q_w;
            //     ::mars::ViewMatrixTextureC<Real, n_qp, dim> q_points = quad.q_p;

            //     auto d_handler = sp.get_dof_handler();

            //     fe_.iterate(MARS_LAMBDA(const Integer elem_index) {
            //         Real gi[dim], gj[dim];
            //         Real J_inv_e[dim * dim];
            //         Real pk[dim];

            //         for (int k = 0; k < dim * dim; ++k) {
            //             J_inv_e[k] = J_inv(elem_index, k);
            //         }

            //         for (int i = 0; i < FEM::elem_nodes; i++) {
            //             const Integer local_dof_i = fe.get_elem_local_dof(elem_index, i);
            //             if (!d_handler.is_owned(local_dof_i)) continue;

            //             for (int j = 0; j < FEM::elem_nodes; j++) {
            //                 const Integer local_dof_j = fe.get_elem_local_dof(elem_index, j);
            //                 // for each dof get the local number

            //                 Real val = 0.0;
            //                 for (int k = 0; k < n_qp; ++k) {
            //                     for (int d = 0; d < dim; ++d) {
            //                         pk[d] = q_points(k, d);
            //                     }

            //                     Elem::Grad::affine_f(i, J_inv_e, pk, gi);
            //                     Elem::Grad::affine_f(j, J_inv_e, pk, gj);

            //                     const Real det_J_e = det_J(elem_index);
            //                     assert(det_J_e > 0.0);
            //                     const Real dx = det_J_e * q_weights(k);

            //                     val += ::mars::Algebra<dim>::dot(gj, gi) * dx;
            //                 }

            //                 sp.atomic_add_value(local_dof_i, local_dof_j, val);
            //             }
            //         }
            //     });
            // }

            // void add_dof_contributions(const ViewMatrixType &res) {
            //     auto dof_handler = sparsity_pattern_.get_dof_handler();
            //     auto eld = fe_.get_elem_dof_map();

            //     dof_handler.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
            //         // update output
            //         for (int i = 0; i < FEM::elem_nodes; i++) {
            //             const Integer local_dof = eld(elem_index, i);
            //             /* atomically updated the contributions to the same dof */
            //             Kokkos::atomic_add(&dof_data(local_dof), res(elem_index, i));
            //         }
            //     });
            // }

            // form the matrix free operator
            // template <Integer INPUT, Integer OUTPUT>
            // void form_operator() {
            //     ViewMatrixType<Real> res("res", sparsity_pattern_.get_dof_handler().get_elem_size(),
            //     FEM::elem_nodes);

            //     integrate<INPUT>(sparsity_pattern_, fe_, quad, det_J_, inv_J_, res);
            //     add_dof_contributions<OUTPUT>(res);
            // }

            // template <class F>
            // void integrate_rhs(F f, ViewMatrixType &res) {
            //     auto det_J = det_J_;
            //     auto sp = sparsity_pattern_.get_dof_handler();
            //     auto fe = fe_;

            //     sp.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
            //         using Elem = typename SparsityPattern::DofHandler::simplex_type;
            //         Real p[2];

            //         const T detj = det_J(elem_index);

            //         for (int i = 0; i < FEM::elem_nodes; i++) {
            //             // forach dof get the local number
            //             const Integer local_dof = fe.get_elem_local_dof(elem_index, i);
            //             sp.template get_dof_coordinates_from_local<Elem::ElemType>(local_dof, p);

            //             const T val = f(p);
            //             const T scaled_val = val * detj / FEM::elem_nodes;
            //             res(elem_index, i) += scaled_val;
            //             /* const Integer owned_index = sp.local_to_owned(local_dof);
            //             Kokkos::atomic_add(&rhs(owned_index), scaled_val); */
            //         }
            //     });
            // }

            // template <class F, Integer RHS>
            // void assemble_local_rhs(F f) {
            //     ViewMatrixType<DMDataType<RHS>> res(
            //         "res", sparsity_pattern_.get_dof_handler().get_elem_size(), FEM::elem_nodes);

            //     integrate_rhs(f, res);
            //     add_dof_contributions<RHS>(res);
            // }

            void init() {
                using Elem = typename SparsityPattern::DofHandler::simplex_type;

                det_J_ = ViewVectorType("detJ", fe_.get_fe_dof_map_size());
                inv_J_ = ViewMatrixType("J_inv", fe_.get_fe_dof_map_size(), Dim * Dim);
                compute_invJ_and_detJ(sparsity_pattern_, fe_, det_J_, inv_J_);

                quad = Quadrature::make();
            }

            MARS_INLINE_FUNCTION SparsityPattern &sp() { return sparsity_pattern_; }
            MARS_INLINE_FUNCTION const SparsityPattern &sp() const { return sparsity_pattern_; }

            MARS_INLINE_FUNCTION ViewVectorType det_J() const { return det_J_; }
            MARS_INLINE_FUNCTION ViewMatrixType J_inv() const { return inv_J_; }
            /* MARS_INLINE_FUNCTION Quadrature quad() const { return
             * quad; } */

        private:
            SparsityPattern &sparsity_pattern_;
            FEM &fe_;

            ViewVectorType det_J_;
            ViewMatrixType inv_J_;
            Quadrature quad;
        };

    }  // namespace mars
}  // namespace utopia
#endif

#ifndef UTOPIA_MARS_FE_HPP
#define UTOPIA_MARS_FE_HPP

#include "mars.hpp"
namespace utopia {
    namespace mars {

        /* using SparsityPattern = DM<DistributedQuad4Mesh, 1, Real, Real>; */
        // template <typename T, typename Scalar>
        // MARS_INLINE_FUNCTION bool invert3(const T *mat, T *mat_inv, const Scalar
        // &det) {
        //   assert(det != 0.);

        //   if (det == 0.) {
        //     return false;
        //   }

        //   mat_inv[0] = (mat[4] * mat[8] - mat[5] * mat[7]) / det;
        //   mat_inv[1] = (mat[2] * mat[7] - mat[1] * mat[8]) / det;
        //   mat_inv[2] = (mat[1] * mat[5] - mat[2] * mat[4]) / det;
        //   mat_inv[3] = (mat[5] * mat[6] - mat[3] * mat[8]) / det;
        //   mat_inv[4] = (mat[0] * mat[8] - mat[2] * mat[6]) / det;
        //   mat_inv[5] = (mat[2] * mat[3] - mat[0] * mat[5]) / det;
        //   mat_inv[6] = (mat[3] * mat[7] - mat[4] * mat[6]) / det;
        //   mat_inv[7] = (mat[1] * mat[6] - mat[0] * mat[7]) / det;
        //   mat_inv[8] = (mat[0] * mat[4] - mat[1] * mat[3]) / det;
        //   return true;
        // }

        // q.points = {{0.5, 0.5},
        // {0.98304589153964795245728880523899, 0.5},
        // {0.72780186391809642112479237299488, 0.074042673347699754349082179816666},
        // {0.72780186391809642112479237299488, 0.92595732665230024565091782018333},
        // {0.13418502421343273531598225407969, 0.18454360551162298687829339850317},
        // {0.13418502421343273531598225407969, 0.81545639448837701312170660149683}};

        // q.weights = {0.28571428571428571428571428571428,
        // 0.10989010989010989010989010989011,
        // 0.14151805175188302631601261486295,
        // 0.14151805175188302631601261486295,
        // 0.16067975044591917148618518733485,
        // 0.16067975044591917148618518733485};

        template <typename T, typename Scalar>
        MARS_INLINE_FUNCTION bool invert2(const T *mat, T *mat_inv, const Scalar &det) {
            mat_inv[0] = mat[3] / det;
            mat_inv[1] = -mat[1] / det;
            mat_inv[2] = -mat[2] / det;
            mat_inv[3] = mat[0] / det;
            return true;
        }

        template <class SparsityPattern, class FEM>
        class FE {
        public:
            template <Integer idx>

            // use as more readable tuple index to identify the data
            FE(SparsityPattern &d, FEM &f) : sparsity_pattern_(d), fe_(f) {}

            // OK as long as we are doing affine meshes in 2D
            template <Integer Type>
            void compute_invJ_and_detJ(const SparsityPattern &data,
                                       const FEM &fe,
                                       ViewVectorType<Real> detJ,
                                       ViewMatrixType<Real> invJ) {
                auto dm = data.get_dof_handler();
                dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
                    Real J[4];
                    Real J_inv[4];
                    Real point_ref[2];
                    Real next_point[2];

                    Integer local_dof = fe.get_elem_local_dof(elem_index, 0);
                    dm.template get_dof_coordinates_from_local<Type>(local_dof, point_ref);

                    local_dof = fe.get_elem_local_dof(elem_index, 1);
                    dm.template get_dof_coordinates_from_local<Type>(local_dof, next_point);

                    // col 0, p1
                    J[0] = next_point[0] - point_ref[0];
                    J[2] = next_point[1] - point_ref[1];

                    // we skip p2 and get p3
                    local_dof = fe.get_elem_local_dof(elem_index, 3);
                    dm.template get_dof_coordinates_from_local<Type>(local_dof, next_point);

                    // col 1, p3
                    J[1] = next_point[0] - point_ref[0];
                    J[3] = next_point[1] - point_ref[1];

                    // compute the determinant and the inverse of the Jacobian
                    Real det_J;
                    Invert<2>::apply(J, J_inv, det_J);
                    detJ(elem_index) = Kokkos::ArithTraits<Real>::abs(det_J);
                    for (int i = 0; i < 4; i++) {
                        invJ(elem_index, i) = J_inv[i];
                    }
                });

                Real measure = KokkosBlas::nrm1(det_J_);
                std::cout << "measure: " << measure << std::endl;
            }

            template <Integer INPUT>
            MARS_INLINE_FUNCTION void gather_elem_data(const SparsityPattern &dm,
                                                       const FEM &fe,
                                                       const Integer elem_index,
                                                       DMDataType<INPUT> *sol) {
                for (int i = 0; i < FEM::elem_nodes; i++) {
                    // forach dof get the local number
                    const Integer local_dof = fe.get_elem_local_dof(elem_index, i);
                    // use the local number to read the corresponding user data
                    sol[i] = dm.template get_dof_data<INPUT>(local_dof);
                }
            }

            // if non-linear than the quad rule should be computed per quad point and not
            // anymore for each element so the coalescing of the Jinv will not matter.
            // In that case maybe a better way to go is parallel through the quad points.
            template <Integer INPUT>
            void integrate(const SparsityPattern &dm,
                           const FEM &fe,
                           const FEQuad4<Real>::Quadrature &quad,
                           ViewVectorType<Real> det_J,
                           ViewMatrixType<Real> J_inv,
                           ViewMatrixType<Real> res) {
                constexpr int n_qp = FEQuad4<Real>::Quadrature::n_points();
                constexpr int dim = FEQuad4<Real>::Quadrature::dim();

                ViewVectorTextureC<Real, n_qp> q_weights = quad.q_w;
                ViewMatrixTextureC<Real, n_qp, dim> q_points = quad.q_p;

                dm.get_dof_handler().elem_iterate(MARS_LAMBDA(const Integer elem_index) {
                    Real gi[dim], g[dim];
                    Real J_inv_e[dim * dim];
                    Real pk[dim];
                    Real sol[FEM::elem_nodes];

                    for (int k = 0; k < dim * dim; ++k) {
                        J_inv_e[k] = J_inv(elem_index, k);
                    }

                    gather_elem_data<INPUT>(dm, fe, elem_index, sol);

                    for (int k = 0; k < n_qp; ++k) {
                        for (int d = 0; d < dim; ++d) {
                            pk[d] = q_points(k, d);
                        }

                        ////////////////////////
                        // Compute physical gradient of solution once per quadrature point
                        FEQuad4<Real>::Grad::ref(pk, sol, gi);
                        Algebra<dim>::m_t_v_mult(J_inv_e, gi, g);
                        ////////////////////////

                        assert(det_J(elem_index) > 0.0);
                        const Real dx = det_J(elem_index) * q_weights(k);

                        for (int i = 0; i < FEM::elem_nodes; i++) {
                            // for each dof get the local number
                            const Integer local_dof = fe.get_elem_local_dof(elem_index, i);
                            FEQuad4<Real>::Grad::affine_f(i, J_inv_e, pk, gi);

                            res(elem_index, i) += Algebra<dim>::dot(g, gi) * dx;

                            assert(res(elem_index, i) == res(elem_index, i));
                        }
                    }
                });
            }

            template <Integer OUTPUT>
            void add_dof_contributions(const ViewMatrixType<Real> &res) {
                auto dof_handler = sparsity_pattern_.get_dof_handler();
                auto eld = fe_.get_elem_dof_map();
                auto dof_data = sparsity_pattern_.template get_dof_data<OUTPUT>();
                dof_handler.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
                    // update output
                    for (int i = 0; i < FEM::elem_nodes; i++) {
                        const Integer local_dof = eld(elem_index, i);
                        /* atomically updated the contributions to the same dof */
                        Kokkos::atomic_add(&dof_data(local_dof), res(elem_index, i));
                    }
                });
            }

            // form the matrix free operator
            template <Integer INPUT, Integer OUTPUT>
            void form_operator() {
                ViewMatrixType<Real> res("res", sparsity_pattern_.get_dof_handler().get_elem_size(), FEM::elem_nodes);

                integrate<INPUT>(sparsity_pattern_, fe_, quad, det_J_, inv_J_, res);
                add_dof_contributions<OUTPUT>(res);
            }

            template <class F, typename T>
            void integrate_rhs(F f, ViewMatrixType<T> &res) {
                auto det_J = det_J_;
                auto dm = sparsity_pattern_.get_dof_handler();
                auto fe = fe_;

                dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
                    using Elem = typename SparsityPattern::DofHandler::simplex_type;
                    Real p[2];

                    const T detj = det_J(elem_index);

                    for (int i = 0; i < FEM::elem_nodes; i++) {
                        // forach dof get the local number
                        const Integer local_dof = fe.get_elem_local_dof(elem_index, i);
                        dm.template get_dof_coordinates_from_local<Elem::ElemType>(local_dof, p);

                        const T val = f(p);
                        const T scaled_val = val * detj / FEM::elem_nodes;
                        res(elem_index, i) += scaled_val;
                        /* const Integer owned_index = dm.local_to_owned(local_dof);
                        Kokkos::atomic_add(&rhs(owned_index), scaled_val); */
                    }
                });
            }

            template <class F, Integer RHS>
            void assemble_local_rhs(F f) {
                ViewMatrixType<DMDataType<RHS>> res("res", sparsity_pattern_.get_dof_handler().get_elem_size(), FEM::elem_nodes);

                integrate_rhs(f, res);
                add_dof_contributions<RHS>(res);
            }

            void init() {
                using Elem = typename SparsityPattern::DofHandler::simplex_type;

                det_J_ = ViewVectorType<Real>("detJ", sparsity_pattern_.get_dof_handler().get_elem_size());
                inv_J_ = ViewMatrixType<Real>("J_inv", sparsity_pattern_.get_dof_handler().get_elem_size(), 4);
                compute_invJ_and_detJ<Elem::ElemType>(sparsity_pattern_, fe_, det_J_, inv_J_);

                quad = FEQuad4<Real>::Quadrature::make();
            }

            MARS_INLINE_FUNCTION SparsityPattern &dm() { return sparsity_pattern_; }
            MARS_INLINE_FUNCTION const SparsityPattern &dm() const { return sparsity_pattern_; }

            MARS_INLINE_FUNCTION ViewVectorType<Real> det_J() const { return det_J_; }
            MARS_INLINE_FUNCTION ViewMatrixType<Real> J_inv() const { return inv_J_; }
            /* MARS_INLINE_FUNCTION FEQuad4<Real>::Quadrature quad() const { return
             * quad; } */

        private:
            SparsityPattern &sparsity_pattern_;
            FEM &fe_;

            ViewVectorType<Real> det_J_;
            ViewMatrixType<Real> inv_J_;
            FEQuad4<Real>::Quadrature quad;
        };

    }  // namespace mars
}  // namespace utopia
#endif

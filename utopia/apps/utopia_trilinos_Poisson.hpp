#include "utopia_Base.hpp"
#include "utopia_Poisson3D.hpp"
#include "utopia_Views.hpp"
#include "utopia_TestFunctions.hpp"

#ifdef  WITH_TRILINOS

namespace utopia
{
    template<typename Matrix, typename Vector, int Backend = Traits<Vector>::Backend>
    class Poisson {};

    template<typename Matrix, typename Vector>
    class Poisson<Matrix, Vector, TRILINOS> final: virtual public UnconstrainedExtendedTestFunction<Matrix, Vector>,
                                                    virtual public ConstrainedExtendedTestFunction<Matrix, Vector> 
    {
    public:
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        static const int Dim      = 2;
        static const int NDofs    = 4;
        static const int NQPoints = 6;


        //FIXME
        typedef Kokkos::TeamPolicy<>               TeamPolicy;
        typedef Kokkos::TeamPolicy<>::member_type  MemberType;

        //FIXME
        using DeviceMatrix      = utopia::MatrixView<Kokkos::View<Scalar **>>;
        using DeviceVector      = utopia::VectorView<Kokkos::View<Scalar *>>;

        using DofView           = Kokkos::View<SizeType **>;
        using PointView         = Kokkos::View<Scalar **>;
        using PhysicalQPointView  = Kokkos::View<Scalar ***>;
        using JacobianDeterminantView = Kokkos::View<Scalar *>;
        using JacobianView      = Kokkos::View<Scalar ***>;
        using GradView          = Kokkos::View<Scalar ****>;
        using ElementMatrixView = Kokkos::View<Scalar ***>;

        //reference element
        using RefFunView   = Kokkos::DualView<Scalar **>;
        using RefGradView  = Kokkos::DualView<Scalar ***>;
        using QPointView   = Kokkos::DualView<Scalar **>;
        using QWeightView  = Kokkos::DualView<Scalar *>;
        using Dev = typename Traits<Vector>::Device;


        // using VectorD = utopia::VectorView<Kokkos::View<Scalar[Dim]>>;

        static SizeType compute_n_elements(const SizeType &n)
        {
            SizeType ret = n;
            for(int i = 1; i < Dim; ++i) {
                ret *= n;
            }

            return ret;
        }

        static SizeType compute_n_points(const SizeType &n)
        {
            SizeType ret = n + 1;
            for(int i = 1; i < Dim; ++i) {
                ret *= (n + 1);
            }

            return ret;
        }

        Poisson(const SizeType &n) :
        n_(n),
        n_elements_(compute_n_elements(n)),
        n_points_(compute_n_points(n)),
        dof_("dof", n_elements_, NDofs),
        point_("point", n_points_, Dim),
        physical_q_points_("physical_q_points", n_elements_, NQPoints, Dim),
        jacobian_("jacobian", n_elements_, Dim, Dim),
        jacobian_inverse_("jacobian_inverse", n_elements_, Dim, Dim),
        jacobian_determinant_("jacobian_determinant", n_elements_),
        grad_("grad", n_elements_, NDofs, NQPoints, Dim),
        element_matrix_("mat", n_elements_, NDofs, NDofs),
        //reference element
        ref_fun_("ref_fun",   NDofs, NQPoints),
        ref_grad_("ref_grad", NDofs, NQPoints, Dim, Dim),
        q_points_("q_points", NQPoints, Dim),
        q_weights_("q_weights", NQPoints)
        {
           init();
            // assemble_mass_matrix();
        }

        ~Poisson()
        {}

        void reinit()
        {
            init();
        }

        void get_A_rhs(Matrix &A, Vector &rhs) const
        {

        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            return false;
        }

        bool hessian(const Vector & /*x*/, Matrix &hessian) const override
        {
            // YES, wrap is more effiicient, but we do not want to own matrix ....
            // as RMTR, needs to modify hessian ...
            // wrap(snes_->jacobian, hessian);
            hessian = laplacian_;
            return false;
        }

        bool value(const Vector &x, Scalar &result) const override
        {
            Vector Ax = laplacian_ * x;
            result = 0.5 * dot(Ax, x) + dot(x, rhs_);
            return false;
        }

        Vector initial_guess() const override
        {
            return initial_guess_;
        }

        const Vector &exact_sol() const override
        {
            return exact_sol_;
        }

        Scalar min_function_value() const override
        {
            return -1.013634375000014e+01;
        }

        std::string name() const override
        {
            return "Poisson_Kokkos";
        }

        SizeType dim() const override
        {
            return n_*n_*n_;
        }

        bool exact_sol_known() const override
        {
            return true;
        }

        bool parallel() const override
        {
            return true;
        }

        bool upper_bound(Vector &ub) const
        {
            ub = local_values(n_points_, 0.45);
            return true;
        }

        bool lower_bound(Vector &lb) const
        {
            lb = local_values(n_points_, -9e9);
            return true;
        }

        bool has_upper_bound() const override
        {
            return false;
        }

        bool has_lower_bound() const override
        {
            return false;
        }

        void describe() const override
        {
            std::cout << "n_elements: " << n_elements_ << std::endl;
            std::cout << "n: " << n_ << std::endl;
            disp(rhs_);
        }

        inline const Matrix &laplacian() const
        {
            return laplacian_;
        }

        inline const Vector &rhs() const
        {
            return rhs_;
        }

        // template<typename Scalar, int Dim>
        // class UniformGridElement {
        // public:
        //     ArrayView<SizeType, Dim> idx;
        // };

        // template<typename Scalar, int Dim>
        // class UniformGrid {
        // public:

        //     using DofView   = Kokkos::View<SizeType **>;
        //     using PointView = Kokkos::View<Scalar **>;
        //     using DimView   = Kokkos::DualView<SizeType[Dim]>;
        //     using Point     = utopia::StaticVector<Scalar, Dim>;
        //     using Jacobian  = utopia::StaticMatrix<Scalar, Dim, Dim>;
        //     using Element   = utopia::UniformGridElement<Dim>;

        //     UTOPIA_INLINE_FUNCTION void element(const SizeType &el_index, UniformGridElement &el)
        //     {
        //         SizeType current = el_index;
        //         const SizeType last = Dim - 1;

        //         for(SizeType i = last; i >= 0; --i) {
        //             const SizeType next = current / dims[i];
        //             el.idx[i] = current - next * dims[i];
        //             current = next;
        //         }
        //     }

        //     UTOPIA_INLINE_FUNCTION void point(const SizeType &i, const SizeType &j, Point &p) const
        //     {
        //         for(int d = 0; d < Dim; ++d) {
        //             const Scalar h = 1./size_[d];
        //             p(d) = i * h;
        //         }
        //     }

        //     UTOPIA_INLINE_FUNCTION void jacobian(const SizeType &element_id, Jacobian &J) const
        //     {
        //         J.set(0.0);

        //         for(int d = 0; d < Dim; ++d) {
        //             const Scalar h = 1./size_[d];
        //             J.set(d, d, h);
        //         }
        //     }

        //     UTOPIA_INLINE_FUNCTION void jacobian_inverse(const SizeType &element_id, Jacobian &J) const
        //     {
        //         J.set(0.0);

        //         for(int d = 0; d < Dim; ++d) {
        //             const Scalar h = 1./size_[d];
        //             J.set(d, d, 1./h);
        //         }
        //     }

        //     UTOPIA_INLINE_FUNCTION Scalar jacobian_determinant(const SizeType &element_id) const
        //     {
        //         return J_det_;
        //     }

        //     UniformGrid(const DimView &size) : size_(size)
        //     {
        //         init();
        //     }

        // private:
        //     DimView size_;
        //     Scalar J_det_;

        //     void init()
        //     {
        //         auto host_size = size_.view_host();

        //         J_det_ = 1./host_size[d];
        //         for(int d = 1; d < Dim; ++d) {
        //             const Scalar h = 1./host_size[d];
        //             J_det_ *= h;
        //         }
        //     }
        // };


        // class FunctionSpace {
        // public:

        // };

        // class FiniteElement {
        // public:

        // };

    private:
        SizeType n_;
        SizeType n_elements_;
        SizeType n_points_;
        Vector exact_sol_, initial_guess_, rhs_;

        DofView dof_;
        PointView point_;
        PhysicalQPointView physical_q_points_;
        JacobianView jacobian_, jacobian_inverse_;
        JacobianDeterminantView jacobian_determinant_;
        GradView grad_;


        ElementMatrixView element_matrix_;

        //reference element
        RefFunView ref_fun_;
        RefGradView ref_grad_;
        QPointView  q_points_;
        QWeightView q_weights_;

        Matrix laplacian_, mass_matrix_;

        void init_vectors()
        {
            exact_sol_ = zeros(n_points_);
            initial_guess_ = zeros(n_points_);
            rhs_ = zeros(n_points_);
        }

        void assemble_matrix(Matrix &mat)
        {
            if(empty(mat)) {
                mat = sparse(n_points_, n_points_, NDofs);
                //First time we construct the graph on the cpu

                Write<Matrix> w_(mat, utopia::GLOBAL_ADD);
                for(SizeType k = 0; k < n_elements_; ++k) {
                    auto dof_k  = Kokkos::subview(dof_, k, Kokkos::ALL());
                    auto el_mat = Kokkos::subview(element_matrix_, k, Kokkos::ALL(), Kokkos::ALL());

                    for(SizeType i = 0; i < NDofs; ++i) {
                        for(SizeType j = 0; j < NDofs; ++j) {
                            mat.c_add(dof_k(i), dof_k(j), el_mat(i, j));
                        }
                    }
                }

            } else {
                mat *= 0.;

                auto device_mat = device_view(mat);

                Dev::parallel_for(n_elements_, UTOPIA_LAMBDA(const SizeType &k) {
                    auto dof_k  = Kokkos::subview(dof_, k, Kokkos::ALL());
                    auto el_mat = Kokkos::subview(element_matrix_, k, Kokkos::ALL(), Kokkos::ALL());

                    for(SizeType i = 0; i < NDofs; ++i) {
                        for(SizeType j = 0; j < NDofs; ++j) {
                            device_mat.atomic_add(dof_k(i), dof_k(j), el_mat(i, j));
                        }
                    }
                });
            }

        }

        void init_q_rule()
        {
            auto host_q_points  = q_points_.view_host();
            auto host_q_weights = q_weights_.view_host();
            auto host_fun       = ref_fun_.view_host();
            auto host_grad      = ref_grad_.view_host();

            host_q_points(0, 0) = 0.5;
            host_q_points(0, 1) = 0.5;
            host_q_points(1, 0) = 0.9304589153964795245728880523899,
            host_q_points(1, 1) = 0.5;
            host_q_points(2, 0) = 0.72780186391809642112479237299488;
            host_q_points(2, 1) = 0.074042673347699754349082179816666;
            host_q_points(3, 0) = 0.72780186391809642112479237299488;
            host_q_points(3, 1) = 0.92595732665230024565091782018333;
            host_q_points(4, 0) = 0.13418502421343273531598225407969;
            host_q_points(4, 1) = 0.18454360551162298687829339850317;
            host_q_points(5, 0) = 0.13418502421343273531598225407969;
            host_q_points(5, 1) = 0.81545639448837701312170660149683;



            host_q_weights(0) = 0.28571428571428571428571428571428;
            host_q_weights(1) = 0.10989010989010989010989010989011;
            host_q_weights(2) = 0.14151805175188302631601261486295;
            host_q_weights(3) = 0.14151805175188302631601261486295;
            host_q_weights(4) = 0.16067975044591917148618518733485;
            host_q_weights(5) = 0.16067975044591917148618518733485;

            for(SizeType i = 0; i < NQPoints; ++i) {
                const Scalar x = host_q_points(i, 0);
                const Scalar y = host_q_points(i, 1);

                //fun
                host_fun(0, i) = (1 - x) * (1 - y);
                host_fun(1, i) = x * (1 - y);
                host_fun(2, i) = x*y;
                host_fun(3, i) = (1 - x) * y;

                //grad
                host_grad(0, i, 0) = y - 1.;
                host_grad(0, i, 1) = x - 1.;

                host_grad(1, i, 0) = 1 - y;
                host_grad(1, i, 1) = -x;

                host_grad(2, i, 0) = y;
                host_grad(2, i, 1) = x;

                host_grad(3, i, 0) = -y;
                host_grad(3, i, 1) = (1 - x);
            }

            //FIXME?
            // q_points_.sync<DeviceType>();
            // q_weights_.sync<DeviceType>();
            // ref_fun_.sync<DeviceType>();
            // ref_grad_.sync<DeviceType>();
        }

        void init()
        {
            Chrono c;
            c.start();
            init_mesh();
            init_vectors();
            assemble_laplacian();

            assemble_rhs(UTOPIA_LAMBDA(const DeviceVector &p) {
                return std::sin(20*p.get(0) + 10*p.get(1)) * std::cos(30*p.get(0) + 6*p.get(1));
            });

            init_boundary_conditions();
            c.stop();

            std::cout << "Poisson::init() " << c << std::endl;
        }

        void init_boundary_conditions()
        {
            Vector bc_markers = values(n_points_, 0.0);
            Vector bc_values  = values(n_points_, 0.0);

            auto bcm_view = device_view(bc_markers);
            // auto bcv_view = device_view(bc_values);

            Dev::parallel_for(n_ + 1, UTOPIA_LAMBDA(const SizeType &i)
            {
                //lower boundary
                bcm_view.set(i, 1.0);

                //left boundary
                bcm_view.set(i * (n_ + 1), 1.0);

                //right boundary
                bcm_view.set(i * (n_ + 1) + n_, 1.0);

                //top boundary
                bcm_view.set((n_ + 1) * n_ + i, 1.0);
            });

            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, bc_values);

            const auto &index = this->get_indices_related_to_BC();
            set_zero_rows(laplacian_, index, 1.);

            rhs_ -= e_mul(bc_markers, rhs_);
            rhs_ += bc_values;
        }

        template<class Fun>
        void assemble_rhs(Fun fun)
        {
            auto device_q_points   = q_points_.view_device();
            auto device_fun        = ref_fun_.view_device();
            auto device_qp_weights = q_weights_.view_device();

            Dev::parallel_for(
                n_elements_,
                UTOPIA_LAMBDA(const SizeType &e_id)
                {
                    DeviceVector vec(Kokkos::subview(element_matrix_, e_id, 0, Kokkos::ALL()));
                    vec.set(0.0);
                    const Scalar det_J = jacobian_determinant_(e_id);

                    DeviceMatrix J(Kokkos::subview(jacobian_, e_id, Kokkos::ALL(), Kokkos::ALL()));

                    for(SizeType i = 0; i < NDofs; ++i) {
                        Scalar val = 0.0;

                        for(SizeType k = 0; k < NQPoints; ++k) {
                            const Scalar dx = det_J * device_qp_weights(k);

                            DeviceVector p(Kokkos::subview(physical_q_points_, e_id, k, Kokkos::ALL()));
                            val += fun(p) * device_fun(i, k) * dx;
                        }

                        vec.set(i, val);
                    }
                }
            );

            assemble_vector(rhs_);
        }

        void assemble_vector(Vector &vec)
        {
            if(empty(vec)) {
                vec = zeros(n_points_);
            } else {
                vec *= 0.;
            }

            Write<Vector> w_(vec, utopia::GLOBAL_ADD);
            for(SizeType k = 0; k < n_elements_; ++k) {
                const auto dof_k  = Kokkos::subview(dof_, k, Kokkos::ALL());
                const auto el_vec = Kokkos::subview(element_matrix_, k, 0, Kokkos::ALL());

                for(SizeType i = 0; i < NDofs; ++i) {
                    const auto dof_I = dof_k(i);
                    const auto val_I = el_vec(i);
                    vec.c_add(dof_I, val_I);
                }
            }
        }

        void assemble_mass_matrix()
        {
            auto q_weights_device = q_weights_.view_device();
            auto device_q_points  = q_points_.view_device();
            auto ref_fun_device   = ref_fun_.view_device();

            Dev::parallel_for(
                n_elements_,
                UTOPIA_LAMBDA(const SizeType &e_id)
                {
                    // const SizeType e_id = team_member.league_rank();
                    DeviceMatrix mat(Kokkos::subview(element_matrix_, e_id, Kokkos::ALL(), Kokkos::ALL()));
                    mat.set(0.0);
                    const Scalar det_J = jacobian_determinant_(e_id);

                    for(SizeType i = 0; i < NDofs; ++i) {

                        for(SizeType k = 0; k < NQPoints; ++k) {
                            const Scalar meas = q_weights_device(k) * det_J;
                            const Scalar fi = ref_fun_device(i, k);

                            mat.add(i, i, fi * fi * meas );

                            for(SizeType j = i+1; j < NDofs; ++j) {
                                const Scalar fj = ref_fun_device(j, k);
                                const Scalar v = fi * fj * meas;
                                mat.add(i, j, v);
                                mat.add(j, i, v);
                            }
                        }
                    }

                });

            assemble_matrix(mass_matrix_);
        }

        void assemble_laplacian()
        {
            Chrono c;

            auto q_weights_device = q_weights_.view_device();

            UTOPIA_TRACE_REGION_BEGIN("Poisson::assemble_laplacian");

            c.start();
            Dev::parallel_for(
                n_elements_,
                UTOPIA_LAMBDA(const SizeType &e_id)
                {
                    DeviceMatrix mat(Kokkos::subview(element_matrix_, e_id, Kokkos::ALL(), Kokkos::ALL()));
                    mat.set(0.0);
                    const Scalar det_J = jacobian_determinant_(e_id);

                    for(SizeType i = 0; i < NDofs; ++i) {
                        for(SizeType k = 0; k < NQPoints; ++k) {
                            const Scalar meas = q_weights_device(k) * det_J;

                            DeviceVector g_i(Kokkos::subview(grad_, e_id, i, k, Kokkos::ALL()));
                            mat.add(i, i, dot(g_i, g_i) * meas );

                            for(SizeType j = i+1; j < NDofs; ++j) {
                                DeviceVector g_j(Kokkos::subview(grad_, e_id, j, k, Kokkos::ALL()));

                                const Scalar v = dot(g_i, g_j) * meas;
                                mat.add(i, j, v);
                                mat.add(j, i, v);
                            }
                        }
                    }
                });

            UTOPIA_TRACE_REGION_END("Poisson::assemble_laplacian");
            c.stop();

            std::cout << "Poisson::assemble_laplacian " << c << std::endl;

            assemble_matrix(laplacian_);
        }

        void init_mesh()
        {
            init_q_rule();

            auto device_q_points  = q_points_.view_device();
            auto device_q_weights = q_weights_.view_device();
            auto device_fun       = ref_fun_.view_device();
            auto device_grad      = ref_grad_.view_device();

            const Scalar h = 1./n_;
            Dev::parallel_for(
                n_,
                UTOPIA_LAMBDA(const SizeType &i)
                {
                    // Kokkos::Dev::parallel_for(Kokkos::TeamThreadRange(team_member, n_), [&] (const SizeType j) {
                    for(SizeType j = 0; j < n_; ++j) {

                        const SizeType e_id = i * n_ + j;
                        DeviceVector p(Kokkos::subview(point_, e_id, Kokkos::ALL()));

                        p.set(0, i * h);
                        p.set(1, j * h);

                        const SizeType n_p = n_ + 1;

                        dof_(e_id, 0) = i * n_p + j;
                        dof_(e_id, 1) = i * n_p + (j + 1);
                        //flipped for dof-consistency and ccw local
                        dof_(e_id, 2) = (i + 1) * n_p + (j + 1);
                        dof_(e_id, 3) = (i + 1) * n_p + j;

                        DeviceMatrix J(Kokkos::subview(jacobian_, e_id, Kokkos::ALL(), Kokkos::ALL()));
                        DeviceMatrix J_inv(Kokkos::subview(jacobian_inverse_, e_id, Kokkos::ALL(), Kokkos::ALL()));

                        J.set(0.0);
                        J.set(0, 0, h);
                        J.set(1, 1, h);

                        J_inv.set(0.0);
                        J_inv.set(0, 0, 1./h);
                        J_inv.set(1, 1, 1./h);

                        jacobian_determinant_(e_id) = h*h;

                        for(SizeType k = 0; k < NQPoints; ++k) {
                            DeviceVector qp(Kokkos::subview(device_q_points, k, Kokkos::ALL()));
                            DeviceVector physical_qp(Kokkos::subview(physical_q_points_, e_id, k, Kokkos::ALL()));

                            physical_qp = J * qp;
                            physical_qp += p;
                        }

                        for(SizeType l = 0; l < NDofs; ++l) {
                            for(SizeType k = 0; k < NQPoints; ++k) {
                                DeviceVector g(Kokkos::subview(device_grad, l, k, Kokkos::ALL()));
                                DeviceVector J_inv_g(Kokkos::subview(grad_, e_id, l, k, Kokkos::ALL()));

                                J_inv_g = J_inv * g;
                            }
                        }
                    }//);
                }
            );
        }
    };
}

#endif //WITH_PETSC

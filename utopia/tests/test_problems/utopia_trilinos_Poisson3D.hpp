#include "utopia_Base.hpp"
#include "utopia_Poisson3D.hpp"

#ifdef  WITH_TRILINOS

namespace utopia
{
    template<typename Matrix, typename Vector, int Backend = Traits<Vector>::Backend>
    class Poisson {};

    template<typename Matrix, typename Vector>
    class Poisson<Matrix, Vector, TRILINOS> final :
        virtual public UnconstrainedExtendedTestFunction<Matrix, Vector>,
        virtual public ConstrainedExtendedTestFunction<Matrix, Vector> {

    public:
        static const int Dim      = 2;
        static const int NDofs    = 4;
        static const int NQPoints = 6;

        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        //FIXME
        typedef Kokkos::TeamPolicy<>               TeamPolicy;
        typedef Kokkos::TeamPolicy<>::member_type  MemberType;

        //FIXME
        using DeviceMatrix      = utopia::MatrixView<Kokkos::View<Scalar **>>;
        using DeviceVector      = utopia::VectorView<Kokkos::View<Scalar *>>;

        using DofView           = Kokkos::View<SizeType **>;
        using PointView         = Kokkos::View<Scalar **>;
        using JacobianDeterminantView = Kokkos::View<Scalar *>;
        using JacobianView      = Kokkos::View<Scalar ***>;
        using GradView          = Kokkos::View<Scalar ****>;
        using ElementMatrixView = Kokkos::View<Scalar ***>;

        //reference element
        using RefFunView   = Kokkos::DualView<Scalar **>;
        using RefGradView  = Kokkos::DualView<Scalar ***>;
        using QPointView   = Kokkos::DualView<Scalar **>;
        using QWeightView  = Kokkos::DualView<Scalar *>;

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
            init_mesh();
            assemble_laplace_element_matrices();
            // assemble_mass_element_matrices();

        }

        ~Poisson()
        {}

        void get_A_rhs(Matrix &A, Vector &rhs) const
        {

        }

        bool gradient_no_rhs(const Vector &x, Vector &g) const override
        {

            return false;
        }

        bool hessian(const Vector & /*x*/, Matrix &hessian) const override
        {
            // YES, wrap is more effiicient, but we do not want to own matrix ....
            // as RMTR, needs to modify hessian ...
            // wrap(snes_->jacobian, hessian);

            return false;
        }

        bool value(const Vector &x, Scalar &result) const override
        {

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

        bool upper_bound(Vector &ub) const override
        {
            assert(false);
            return true;
        }

        bool lower_bound(Vector &lb) const override
        {
            assert(false);
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

        void describe() const
        {
            std::cout << "n_elements: " << n_elements_ << std::endl;
            std::cout << "n: " << n_ << std::endl;

            // for(SizeType i = 0; i < n_elements_; ++i) {
            //     DeviceMatrix m(Kokkos::subview(element_matrix_, i, Kokkos::ALL(), Kokkos::ALL()));
            //     disp(m);
            // }

            disp(laplacian_);
        }


    private:
        SizeType n_;
        SizeType n_elements_;
        SizeType n_points_;
        Vector exact_sol_, initial_guess_;

        DofView dof_;
        PointView point_;

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

        void assemble_matrix(Matrix &mat)
        {
            if(empty(mat)) {
                mat = sparse(n_points_, n_points_, NDofs);
            } else {
                mat *= 0.;
            }


            Write<Matrix> w_(mat, utopia::GLOBAL_ADD);
            for(SizeType k = 0; k < n_elements_; ++k) {
                auto dof_k = Kokkos::subview(dof_, k, Kokkos::ALL());
                auto el_mat = Kokkos::subview(element_matrix_, k, Kokkos::ALL(), Kokkos::ALL());

                for(SizeType i = 0; i < NDofs; ++i) {
                    for(SizeType j = 0; j < NDofs; ++j) {
                        mat.c_add(dof_k(i), dof_k(j), el_mat(i, j));
                    }
                }
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
            host_q_points(1, 0) = 0.5; 
            host_q_points(1, 1) = 0.0;
            host_q_points(2, 0) = 0.0; 
            host_q_points(2, 1) = 0.5;
            host_q_points(3, 0) = 1.0/6.0; 
            host_q_points(3, 1) = 1.0/6.0;
            host_q_points(4, 0) = 1.0/6.0; 
            host_q_points(4, 1) = 2.0/3.0;
            host_q_points(5, 0) = 2.0/3.0; 
            host_q_points(5, 1) = 1.0/6.0;

            host_q_weights(0) = 1.0/30.0;
            host_q_weights(1) = 1.0/30.0;
            host_q_weights(2) = 1.0/30.0;
            host_q_weights(3) = 0.3;
            host_q_weights(4) = 0.3;
            host_q_weights(5) = 0.3;

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

            //FIXME synch with device?
        }


        void assemble_mass_element_matrices()
        {
            auto q_weights_device = q_weights_.view_device();
            auto ref_fun_device    = ref_fun_.view_device();

            Kokkos::parallel_for(
                "Poisson::assemble_mass_element_matrices",
                TeamPolicy(n_elements_, Kokkos::AUTO),
                KOKKOS_LAMBDA(const MemberType &team_member) 
                {
                    const SizeType e_id = team_member.league_rank();
                    DeviceMatrix mat(Kokkos::subview(element_matrix_, e_id, Kokkos::ALL(), Kokkos::ALL()));
                    mat.set(0.0);
                    const Scalar det_J = jacobian_determinant_(e_id);

                    Kokkos::parallel_for(
                        Kokkos::TeamThreadRange(team_member, NDofs), [&] (const SizeType i)
                        {
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

                        });

                });

            assemble_matrix(mass_matrix_);
        }


        void assemble_laplace_element_matrices()
        {

            auto q_weights_device = q_weights_.view_device();

            Kokkos::parallel_for(
                "Poisson::assemble_laplace_element_matrices",
                TeamPolicy(n_elements_, Kokkos::AUTO),
                KOKKOS_LAMBDA(const MemberType &team_member) 
                {
                    const SizeType e_id = team_member.league_rank();
                    DeviceMatrix mat(Kokkos::subview(element_matrix_, e_id, Kokkos::ALL(), Kokkos::ALL()));
                    mat.set(0.0);
                    const Scalar det_J = jacobian_determinant_(e_id);

                    Kokkos::parallel_for(
                        Kokkos::TeamThreadRange(team_member, NDofs), [&] (const SizeType i)
                        {
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

                        });
                });


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
            Kokkos::parallel_for(
                "Poisson::init_mesh",
                TeamPolicy(n_, Kokkos::AUTO),
                KOKKOS_LAMBDA(const MemberType &team_member) {
                    const SizeType i = team_member.league_rank();

                    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_), [&] (const SizeType j) {
                        const SizeType e_id = i * n_ + j;
                        DeviceVector p(Kokkos::subview(point_, e_id, Kokkos::ALL()));

                        p.set(0, i * h);
                        p.set(1, j * h);

                        dof_(e_id, 0) = e_id;
                        dof_(e_id, 1) = i * n_ + (j + 1);
                        dof_(e_id, 2) = (i + 1) * n_ + (j + 1);
                        dof_(e_id, 3) = (i + 1) * n_ + j;

                        DeviceMatrix J(Kokkos::subview(jacobian_, e_id, Kokkos::ALL(), Kokkos::ALL()));
                        DeviceMatrix J_inv(Kokkos::subview(jacobian_inverse_, e_id, Kokkos::ALL(), Kokkos::ALL()));
                       
                        J.set(0.0);
                        J.set(0, 0, h);
                        J.set(1, 1, h);

                        J_inv.set(0.0);
                        J_inv.set(0, 0, 1./h);
                        J_inv.set(1, 1, 1./h);

                        jacobian_determinant_(e_id) = h*h;

                        
                        for(SizeType l = 0; l < NDofs; ++l) {
                            for(SizeType k = 0; k < NQPoints; ++k) {
                                DeviceVector g(Kokkos::subview(device_grad, l, k, Kokkos::ALL()));
                                DeviceVector J_inv_g(Kokkos::subview(grad_, e_id, l, k, Kokkos::ALL()));
                                
                                J_inv_g = J_inv * g;
                            }
                        }
                    });
                }
            );
        }
    };
}

#endif //WITH_PETSC

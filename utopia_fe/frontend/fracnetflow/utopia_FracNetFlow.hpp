#ifndef UTOPIA_FRAC_NET_FLOW_HPP
#define UTOPIA_FRAC_NET_FLOW_HPP

#include "utopia_CoupledFEFunction.hpp"
#include "utopia_NCFunctionSpace.hpp"
#include "utopia_NLSolve.hpp"
#include "utopia_QuadraticFEFunction.hpp"

#include "utopia_petsc_CrsView.hpp"

namespace utopia {
    template <class FunctionSpace>
    class StabilizeTransport final : public MatrixTransformer<typename Traits<FunctionSpace>::Matrix>,
                                     public PostProcessor<typename Traits<FunctionSpace>::Vector> {
    public:
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using SizeType = typename Traits<FunctionSpace>::SizeType;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Vector = typename Traits<FunctionSpace>::Vector;
        using IndexArray = typename Traits<Vector>::IndexArray;

        StabilizeTransport(const std::shared_ptr<FunctionSpace> &space) : space_(space) {}

        static void create_matrices(Matrix &in_out, Matrix &diff) {
            Matrix mat_t = transpose(in_out);

            /////////////////////////////////////////////////
            diff = in_out;
            diff.transform([&](const SizeType i, const SizeType j, const Scalar &value) {
                if (i == j) {
                    return 0.0;
                } else {
                    const Scalar value_t = mat_t.get(i, j);
                    Scalar max_val = std::max(value, value_t);

                    if (max_val > 0.0) {
                        max_val *= -1.0;
                        return max_val;
                    } else {
                        return 0.0;
                    }
                }
            });

            Vector diag_elem = -1.0 * sum(diff, 1);
            diff += diag(diag_elem);
            /////////////////////////////////////////////////

            in_out += diff;
        }

        void apply(Matrix &mat) override {
            UTOPIA_TRACE_REGION_BEGIN("StabilizeTransport::apply");

            Matrix mat_diff;
            create_matrices(mat, mat_diff);

            UTOPIA_TRACE_REGION_END("StabilizeTransport::apply");
        }

        void apply_to_matrices(Matrix &mass_matrix, Matrix &op) override {
            UTOPIA_TRACE_SCOPE("StabilizeTransport::apply_to_matrices");

            assert(!mass_matrix.empty());
            assert(!op.empty());

            auto op_diff = std::make_shared<Matrix>();
            create_matrices(op, *op_diff);

            // lump mass matrix and extract diff
            auto mass_diff = std::make_shared<Matrix>();
            *mass_diff = mass_matrix;
            mass_matrix.lump();
            *mass_diff -= mass_matrix;

            initialize(make_ref(op), op_diff, make_ref(mass_matrix), mass_diff);
        }

        void set_time(const std::shared_ptr<SimulationTime<Scalar>> &time) override { dt_ = time->delta(); }

        void initialize(const std::shared_ptr<Matrix> &A_corrected,
                        const std::shared_ptr<Matrix> &A_diff,
                        const std::shared_ptr<Matrix> &M_corrected,
                        const std::shared_ptr<Matrix> &M_diff) {
            A_corrected_ = A_corrected;
            A_diff_ = A_diff;
            M_corrected_ = M_corrected;
            M_diff_ = M_diff;

            inverse_mass_vector_ = diag(*M_corrected);
            inverse_mass_vector_ = 1. / inverse_mass_vector_;
        }

        static void add_contrib(Matrix &A, const Vector &x, const Vector &u_ghosts, const Scalar factor, Matrix &out) {
            auto x_d = local_view_host(x);
            auto x_o = local_view_host(u_ghosts);

            PetscCrsView A_d, A_o;
            views_host(A, A_d, A_o);

            PetscCrsView out_d, out_o;
            views_host(out, out_d, out_o);

            SizeType rows = A_d.rows();
            for (SizeType r = 0; r < rows; r++) {
                Scalar xi = x_d.get(r);

                {
                    auto row = A_d.row(r);
                    auto out_row = out_d.row(r);

                    for (SizeType k = 0; k < row.length; k++) {
                        SizeType c = row.colidx(k);
                        Scalar Aij = row.value(k);
                        Scalar xj = x_d.get(c);

                        assert(out_row.colidx(k) == c);
                        out_row.value(k) += factor * Aij * (xi - xj);
                    }
                }

                if (!empty(A_o)) {
                    auto row = A_o.row(r);
                    auto out_row = out_o.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        SizeType c = row.colidx(k);
                        Scalar Aij = row.value(k);
                        Scalar xj = x_o.get(c);

                        assert(out_row.colidx(k) == c);
                        out_row.value(k) += factor * Aij * (xi - xj);
                    }
                }
            }
        }

        // Eq (79)
        static void pre_limiting_step(const Vector &x, const Vector &u_ghosts, Matrix &F) {
            auto x_d = local_view_host(x);
            auto x_o = local_view_host(u_ghosts);

            // const Scalar tol = -1e-8;
            // const Scalar tol = -1e-10;
            const Scalar tol = -1e-12;
            // const Scalar tol = 0;

            PetscCrsView F_d, F_o;
            views_host(F, F_d, F_o);

            SizeType rows = F_d.rows();
            for (SizeType r = 0; r < rows; r++) {
                Scalar xi = x_d.get(r);

                {
                    auto F_row = F_d.row(r);
                    for (SizeType k = 0; k < F_row.length; k++) {
                        SizeType c = F_row.colidx(k);
                        Scalar Fij = F_row.value(k);
                        Scalar xj = x_d.get(c);

                        if (Fij * (xj - xi) > tol) {
                            F_row.value(k) = 0;
                        }
                    }
                }

                if (!empty(F_o)) {
                    auto F_row = F_o.row(r);
                    for (SizeType k = 0; k < F_row.length; k++) {
                        SizeType c = F_row.colidx(k);
                        Scalar Fij = F_row.value(k);
                        Scalar xj = x_o.get(c);

                        if (Fij * (xj - xi) > tol) {
                            F_row.value(k) = 0;
                        }
                    }
                }
            }
        }

        // Book_Kuzmin.pdf page 163 eq (80)
        static void negative_postive_antidiffusive_fluxes(Matrix &F, Vector &P_minus, Vector &P_plus) {
            auto P_plus_d = local_view_host(P_plus);
            auto P_minus_d = local_view_host(P_minus);

            PetscCrsView F_d, F_o;
            views_host(F, F_d, F_o);

            SizeType rows = F_d.rows();
            for (SizeType r = 0; r < rows; r++) {
                Scalar P_plus_i = 0;
                Scalar P_minus_i = 0;

                {
                    auto row = F_d.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        Scalar Fij = row.value(k);
                        P_plus_i += std::max(0., Fij);
                        P_minus_i += std::min(0., Fij);
                    }
                }

                if (!empty(F_o)) {
                    auto row = F_o.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        Scalar Fij = row.value(k);
                        P_plus_i += std::max(0., Fij);
                        P_minus_i += std::min(0., Fij);
                    }
                }

                P_plus_d.set(r, P_plus_i);
                P_minus_d.set(r, P_minus_i);
            }
        }

        // Book_Kuzmin.pdf page 163 eq (81)
        static void min_max_bound(const Scalar dt,
                                  Matrix &M_corrected,
                                  Vector &u,
                                  Vector &u_ghosts,
                                  Vector &Q_minus,
                                  Vector &Q_plus) {
            auto u_d = local_view_host(u);
            auto u_o = local_view_host(u_ghosts);

            auto Q_plus_d = local_view_host(Q_plus);
            auto Q_minus_d = local_view_host(Q_minus);

            PetscCrsView M_d, M_o;
            views_host(M_corrected, M_d, M_o);

            const SizeType rows = M_d.rows();
            for (SizeType r = 0; r < rows; r++) {
                Scalar u_max_i = u_d.get(r);
                Scalar u_min_i = u_d.get(r);
                Scalar Mii = 0;

                {
                    auto row = M_d.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        const SizeType c = row.colidx(k);
                        const Scalar Mij = row.value(k);

                        if (c == r) {
                            Mii = Mij;
                        }

                        u_max_i = std::max(u_max_i, u_d.get(c));
                        u_min_i = std::min(u_min_i, u_d.get(c));
                    }
                }

                if (!empty(M_o)) {
                    auto row = M_o.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        const SizeType c = row.colidx(k);
                        u_max_i = std::max(u_max_i, u_o.get(c));
                        u_min_i = std::min(u_min_i, u_o.get(c));
                    }
                }

                Q_plus_d.set(r, (Mii / dt) * (u_max_i - u_d.get(r)));
                Q_minus_d.set(r, (Mii / dt) * (u_min_i - u_d.get(r)));
            }
        }

        // Eq (66) + (82)
        // check Eq (115)
        // static void create_update(Matrix &F,
        //                           const Vector &R_minus,
        //                           const Vector &R_minus_ghosts,
        //                           const Vector &R_plus,
        //                           const Vector &R_plus_ghosts,
        //                           Vector &f_bar) {
        //     auto R_plus_d = local_view_host(R_plus);
        //     auto R_minus_d = local_view_host(R_minus);

        //     auto R_plus_o = local_view_host(R_plus_ghosts);
        //     auto R_minus_o = local_view_host(R_minus_ghosts);

        //     auto f_bar_o = local_view_host(f_bar);

        //     PetscCrsView F_d, F_o;
        //     views_host(F, F_d, F_o);

        //     SizeType rows = F_d.rows();
        //     for (SizeType r = 0; r < rows; r++) {
        //         Scalar fi = 0;

        //         {
        //             auto row = F_d.row(r);
        //             for (SizeType k = 0; k < row.length; k++) {
        //                 const SizeType c = row.colidx(k);
        //                 if (c == r) continue;

        //                 Scalar Fij = row.value(k);
        //                 Scalar alpha_ij = 0;

        //                 if (Fij > 0) {
        //                     alpha_ij = std::min(R_plus_d.get(r), R_minus_d.get(c));
        //                 } else {
        //                     alpha_ij = std::min(R_minus_d.get(r), R_plus_d.get(c));
        //                 }

        //                 fi += alpha_ij * Fij;
        //             }
        //         }

        //         if (!empty(F_o)) {
        //             auto row = F_o.row(r);
        //             for (SizeType k = 0; k < row.length; k++) {
        //                 const SizeType c = row.colidx(k);

        //                 Scalar Fij = row.value(k);

        //                 Scalar alpha_ij = 0;

        //                 if (Fij > 0) {
        //                     alpha_ij = std::min(R_plus_d.get(r), R_minus_o.get(c));
        //                 } else {
        //                     alpha_ij = std::min(R_minus_d.get(r), R_plus_o.get(c));
        //                 }

        //                 fi += alpha_ij * Fij;
        //             }
        //         }

        //         f_bar_o.set(r, fi);
        //     }
        // }

        void create_update(Matrix &F,
                           const Vector &R_minus,
                           const Vector &R_minus_ghosts,
                           const Vector &R_plus,
                           const Vector &R_plus_ghosts,
                           Vector &f_bar) const {
            auto R_plus_d = local_view_host(R_plus);
            auto R_minus_d = local_view_host(R_minus);

            auto R_plus_o = local_view_host(R_plus_ghosts);
            auto R_minus_o = local_view_host(R_minus_ghosts);

            auto f_bar_o = local_view_host(f_bar);

            PetscCrsView F_d, F_o;
            views_host(F, F_d, F_o);

            Matrix alpha = F;
            alpha *= 0.;

            PetscCrsView alpha_d, alpha_o;
            views_host(alpha, alpha_d, alpha_o);

            SizeType rows = F_d.rows();
            for (SizeType r = 0; r < rows; r++) {
                Scalar fi = 0;

                {
                    auto row = F_d.row(r);
                    auto alpha_row = alpha_d.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        const SizeType c = row.colidx(k);
                        if (c == r) continue;

                        Scalar Fij = row.value(k);
                        Scalar alpha_ij = 0;

                        if (Fij > 0) {
                            alpha_ij = std::min(R_plus_d.get(r), R_minus_d.get(c));
                        } else {
                            alpha_ij = std::min(R_minus_d.get(r), R_plus_d.get(c));
                        }

                        fi += alpha_ij * Fij;

                        alpha_row.value(k) = alpha_ij;
                    }
                }

                if (!empty(F_o)) {
                    auto row = F_o.row(r);
                    auto alpha_row = alpha_o.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        const SizeType c = row.colidx(k);

                        Scalar Fij = row.value(k);

                        Scalar alpha_ij = 0;

                        if (Fij > 0) {
                            alpha_ij = std::min(R_plus_d.get(r), R_minus_o.get(c));
                        } else {
                            alpha_ij = std::min(R_minus_d.get(r), R_plus_o.get(c));
                        }

                        fi += alpha_ij * Fij;
                        alpha_row.value(k) = alpha_ij;
                    }
                }

                f_bar_o.set(r, fi);
            }

            if (debug_) {
                rename("alpha", alpha);
                write("load_alpha.m", alpha);
            }
        }

        void post_process(const Vector &g, Vector &u) const override {
            UTOPIA_TRACE_SCOPE("StabilizeTransport::post_process");
            //!!! CHECK -g instead of g because the rhs is = -g
            Vector u_dot = (*A_corrected_) * u + g;

            // TODO check if +
            // u_dot = -u_dot;

            space_->apply_zero_constraints(u_dot);

            // Vector u_dot = (*A_corrected_) * x + g;
            u_dot = e_mul(inverse_mass_vector_, u_dot);

            if (debug_) {
                space_->write("u_dot.e", u_dot);
                Vector force = e_mul(inverse_mass_vector_, g);
                space_->write("force.e", force);
            }

            Matrix F = *A_corrected_;
            F *= 0;

            IndexArray ghosts;
            A_corrected_->ghosts(ghosts);

            // Split vectors into local block and off diagonal block
            // retrieve ghost values of vectors
            Vector u_ghosts, u_dot_ghosts;
            u.select(ghosts, u_ghosts);
            u_dot.select(ghosts, u_dot_ghosts);

            // TODO check if -1 or 1
            add_contrib(*A_diff_, u, u_ghosts, -1, F);
            add_contrib(*M_diff_, u_dot, u_dot_ghosts, 1, F);
            pre_limiting_step(u, u_ghosts, F);

            if (debug_) {
                Matrix zeros = F + transpose(F);
                Scalar norm_zeros = norm2(zeros);

                if (!u.comm().rank()) {
                    utopia::out() << "norm_zeros: " << norm_zeros << "\n";
                }
            }

            if (debug_) {
                rename("F", F);
                write("load_F.m", F);
            }

            Vector P_plus(layout(u), 0), P_minus(layout(u), 0);
            negative_postive_antidiffusive_fluxes(F, P_plus, P_minus);

            if (debug_) {
                space_->write("P_plus.e", P_plus);
                space_->write("P_minus.e", P_minus);
            }

            Vector Q_plus(layout(u), 0), Q_minus(layout(u), 0);
            min_max_bound(dt_, *M_corrected_, u, u_ghosts, Q_minus, Q_plus);

            if (debug_) {
                space_->write("Q_plus.e", Q_plus);
                space_->write("Q_minus.e", Q_minus);
            }

            // Book_Kuzmin.pdf page 163 eq (82)
            Vector R_plus(layout(u), 0), R_minus(layout(u), 0);
            R_plus = Q_plus / P_plus;
            R_plus.e_min(1);

            R_minus = Q_minus / P_minus;
            R_minus.e_min(1);

            // Dirichlet nodes should not change
            space_->apply_value_constraints(1, R_plus);
            space_->apply_value_constraints(1, R_minus);

            if (debug_) {
                space_->write("R_plus.e", R_plus);
                space_->write("R_minus.e", R_minus);
            }

            Vector R_plus_ghosts, R_minus_ghosts;
            R_plus.select(ghosts, R_plus_ghosts);
            R_minus.select(ghosts, R_minus_ghosts);

            Vector u_correction(layout(u), 0);
            create_update(F, R_minus, R_minus_ghosts, R_plus, R_minus_ghosts, u_correction);

            u_correction *= dt_;
            u_correction = e_mul(inverse_mass_vector_, u_correction);

            space_->apply_zero_constraints(u_correction);

            Scalar norm_u_correction = norm2(u_correction);
            if (!u.comm().rank()) {
                utopia::out() << "norm_u_correction: " << norm_u_correction << "\n";
            }

            u += u_correction;

            if (debug_) {
                space_->write("u_correction.e", u_correction);
                Utopia::Abort("Exiting for dbg purposes!");
            }
        }

    private:
        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<Matrix> A_corrected_;
        std::shared_ptr<Matrix> A_diff_;
        std::shared_ptr<Matrix> M_corrected_;
        std::shared_ptr<Matrix> M_diff_;
        Vector inverse_mass_vector_;
        Scalar dt_{1};
        bool debug_{true};
        // bool debug_{false};
    };

    template <class FunctionSpace>
    class FracNetFlow : public NLSolve<FunctionSpace> {
    public:
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;

        using Super = utopia::NLSolve<FunctionSpace>;
        using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using ImplicitEulerIntegrator_t = utopia::ImplicitEulerIntegrator<FunctionSpace>;
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;
        using TimeDependentFunction_t = utopia::TimeDependentFunction<FunctionSpace>;
        using CoupledFEFunction_t = utopia::CoupledFEFunction<FunctionSpace>;
        using QuadraticFEFunction_t = utopia::QuadraticFEFunction<FunctionSpace>;

        std::shared_ptr<FEFunctionInterface_t> create_problem(
            Input &in,
            const std::string &problem_type,
            const std::shared_ptr<FunctionSpace> &porous_matrix,
            const std::vector<std::shared_ptr<FunctionSpace>> &fracture_networks) {
            auto problem = std::make_shared<CoupledFEFunction_t>();

            std::shared_ptr<FEFunctionInterface_t> matrix_problem;

            if (fracture_networks.empty()) {
                this->warning("fracture_networks is undefined");
            }

            in.get("porous_matrix", [&](Input &node) {
                this->status("Reading material for " + porous_matrix->name());

                auto flow = std::make_shared<FEModelFunction_t>(porous_matrix);
                flow->set_environment(this->environment());

                node.require(problem_type, *flow);
                matrix_problem = flow;
                problem->add_master_function(porous_matrix->name(), matrix_problem);
            });

            in.get("fracture_networks", [&](Input &array_node) {
                int idx = 0;
                array_node.get_all([&](Input &node) {
                    auto space = fracture_networks[idx++];

                    this->status("Reading material for " + space->name());

                    auto flow = std::make_shared<FEModelFunction_t>(space);
                    flow->set_environment(this->environment());
                    node.require(problem_type, *flow);
                    problem->add_function(space->name(), flow);
                    auto &c = problem->add_coupling(porous_matrix->name(), space->name());
                    node.get("coupling", c);
                });
            });

            this->status("Initializing coupled problem");
            problem->initialize();

            std::string integrator;
            in.get("integrator", integrator);

            if (problem_type == "transport") {
                auto st = std::make_shared<StabilizeTransport<FunctionSpace>>(porous_matrix);
                problem->add_matrix_transformer(st);

                // FIXME uncomment after fix!
                problem->add_post_processor(st);
            }

            if (problem_type == "transport" || integrator == "ImplicitEuler") {
                if (simplify_problem_ && problem->is_linear()) {
                    return std::make_shared<ImplicitEulerIntegrator_t>(
                        std::make_shared<QuadraticFEFunction_t>(problem));
                } else {
                    return std::make_shared<ImplicitEulerIntegrator_t>(problem);
                }
            } else {
                return problem;
            }
        }

        std::string name() const override { return "FracNetFlow"; }

        void read(Input &in) override {
            Super::read(in);

            in.get("simplify_problem", simplify_problem_);

            std::shared_ptr<FunctionSpace> porous_matrix;
            std::vector<std::shared_ptr<FunctionSpace>> fracture_networks;

            in.get("porous_matrix", [this, &porous_matrix](Input &node) {
                // Read the function-space of the porous-matrix

                node.get("space", [this, &porous_matrix](Input &space_node) {
                    // auto s = std::make_shared<FunctionSpace>(this->comm());
                    auto s = std::make_shared<NCFunctionSpace<FunctionSpace>>(this->comm());
                    // Use this so everyhting is added to the env automatically when calling read
                    // s->set_environment(this->environment());

                    bool read_state = false;
                    space_node.get("read_state", read_state);
                    if (read_state) {
                        auto field = std::make_shared<Field<FunctionSpace>>();
                        s->read_with_state(space_node, *field);
                        this->environment()->add_field(field);
                    } else {
                        s->read(space_node);
                    }

                    if (s->name().empty()) {
                        utopia::err() << "Value for key \"name\" must be defined for space node\n";
                        Utopia::Abort();
                    }

                    this->environment()->add_space(s);
                    porous_matrix = s;
                });
            });

            if (!porous_matrix) {
                Utopia::Abort("Definition of \"space\" undefined for node porous_matrix!");
            }

            in.get("fracture_networks", [this, &fracture_networks](Input &array_node) {
                array_node.get_all([this, &fracture_networks](Input &node) {
                    // Read the function-space of the fracture-network

                    node.get("space", [this, &fracture_networks](Input &space_node) {
                        // auto s = std::make_shared<FunctionSpace>(this->comm());
                        auto s = std::make_shared<NCFunctionSpace<FunctionSpace>>(this->comm());
                        // Use this so everyhting is added to the env automatically when calling read
                        // s->set_environment(this->environment());

                        bool read_state = false;
                        space_node.get("read_state", read_state);
                        if (read_state) {
                            auto field = std::make_shared<Field<FunctionSpace>>();
                            s->read_with_state(space_node, *field);
                            this->environment()->add_field(field);
                        } else {
                            s->read(space_node);
                        }

                        if (s->name().empty()) {
                            utopia::err() << "Value for key \"name\" must be defined for space node\n";
                            Utopia::Abort();
                        }

                        this->environment()->add_space(s);
                        fracture_networks.push_back(s);
                    });
                });
            });

            std::string problem_type = "flow";
            in.get("problem_type", problem_type);

            auto problem = create_problem(in, problem_type, porous_matrix, fracture_networks);

            problem->read(in);

            this->init(problem);
        }

        bool simplify_problem_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_FRAC_NET_FLOW_HPP

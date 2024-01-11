#ifndef UTOPIA_MATRIX_TRANSFORMER_HPP
#define UTOPIA_MATRIX_TRANSFORMER_HPP

#include <memory>
#include "utopia_Input.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Traits.hpp"
#include "utopia_SimulationTime.hpp"

// #include "utopia_StabilizeTransport.hpp"

namespace utopia {
    template <class Matrix>
    class MatrixTransformer : public Configurable {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;
        using MatrixTransformerPtr = std::unique_ptr<MatrixTransformer>;
        virtual ~MatrixTransformer() = default;
        virtual void apply(Matrix &mat) = 0;
        virtual void apply_to_matrices(Matrix &mass_matrix, Matrix &op) = 0;
        virtual void set_time(const std::shared_ptr<SimulationTime<Scalar>> &time) = 0;
        void read(Input &) override {}

        class Registry {
        public:
            static Registry &instance() {
                static Registry instance_;
                return instance_;
            }

            template <class Type>
            void register_transformer(const std::string &name) {
                transformers[name] = []() -> MatrixTransformerPtr { return utopia::make_unique<Type>(); };
            }

            MatrixTransformerPtr find_transformer(const std::string &name) const {
                MatrixTransformerPtr ret;
                auto it = transformers.find(name);

                if (it != transformers.end()) {
                    ret = it->second();
                }

                return ret;
            }

        private:
            std::map<std::string, std::function<MatrixTransformerPtr()>> transformers;
        };

        static std::unique_ptr<MatrixTransformer> New(const std::string &type) {
            return Registry::instance().find_transformer(type);
        }

        template <class SubType>
        static void Register(const std::string &type) {
            return Registry::instance().template register_transformer<SubType>(type);
        }
    };

    template <class Vector>
    class PostProcessor {
    public:
        virtual ~PostProcessor() = default;
        virtual void post_process(const Vector &rhs, Vector &x) const = 0;
    };

    template <class Matrix>
    class StabilizeTransport final : public MatrixTransformer<Matrix>,
                                     public PostProcessor<typename Traits<Matrix>::Vector> {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;
        using Vector = typename Traits<Matrix>::Vector;
        using IndexArray = typename Traits<Vector>::IndexArray;

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

        void set_time(const std::shared_ptr<SimulationTime<Scalar>> &time) override
        {
            dt_ = time->delta();
        }

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

        static void add_contrib(Matrix &A, const Vector &x, const Vector &x_ghosts, Matrix &out) {
            auto x_view_d = local_view_host(x);
            auto x_view_o = local_view_host(x_ghosts);

            PetscCrsView A_d, A_o;
            views_host(A, A_d, A_o);

            PetscCrsView out_d, out_o;
            views_host(out, out_d, out_o);

            SizeType rows = A_d.rows();
            for (SizeType r = 0; r < rows; r++) {
                Scalar xi = x_view_d.get(r);

                {
                    auto row = A_d.row(r);
                    auto out_row = out_d.row(r);

                    for (SizeType k = 0; k < row.length; k++) {
                        SizeType c = row.colidx(k);
                        Scalar Aij = row.value(k);
                        Scalar xj = x_view_d.get(c);

                        assert(out_row.colidx(k) == c);
                        out_row.value(k) += Aij * (xi - xj);
                    }
                }

                if(!empty(A_o))
                {
                    auto row = A_o.row(r);
                    auto out_row = out_o.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        SizeType c = row.colidx(k);
                        Scalar Aij = row.value(k);
                        Scalar xj = x_view_o.get(c);

                        assert(out_row.colidx(k) == c);
                        out_row.value(k) += Aij * (xi - xj);
                    }
                }
            }
        }

        // Book_Kuzmin.pdf page 163 eq (80)
        static void negative_postive_antidiffusive_fluxes(Matrix &F, Vector &P_minus, Vector &P_plus) {
            auto P_plus_view = local_view_host(P_plus);
            auto P_minus_view = local_view_host(P_minus);

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

                if(!empty(F_o))
                {
                    auto row = F_o.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        Scalar Fij = row.value(k);
                        P_plus_i += std::max(0., Fij);
                        P_minus_i += std::min(0., Fij);
                    }
                }

                P_plus_view.set(r, P_plus_i);
                P_minus_view.set(r, P_minus_i);
            }
        }

        // Book_Kuzmin.pdf page 163 eq (81)
        static void min_max_bound(const Scalar dt,
                                  Matrix &M_corrected,
                                  Vector &u,
                                  Vector &u_ghosts,
                                  Vector &Q_minus,
                                  Vector &Q_plus) {
            auto u_view_d = local_view_host(u);
            auto u_view_o = local_view_host(u_ghosts);

            auto Q_plus_view = local_view_host(Q_plus);
            auto Q_minus_view = local_view_host(Q_minus);

            PetscCrsView M_d, M_o;
            views_host(M_corrected, M_d, M_o);

            const SizeType rows = M_d.rows();
            for (SizeType r = 0; r < rows; r++) {
                Scalar u_max_i = u_view_d.get(r);
                Scalar u_min_i = u_view_d.get(r);
                Scalar Mii = 0;

                {
                    auto row = M_d.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        const SizeType c = row.colidx(k);
                        const Scalar Mij = row.value(k);

                        if (c == r) {
                            Mii = Mij;
                        }

                        u_max_i = std::max(u_max_i, u_view_d.get(c));
                        u_min_i = std::min(u_min_i, u_view_d.get(c));
                    }
                }

                if(!empty(M_o))
                {
                    auto row = M_o.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        const SizeType c = row.colidx(k);
                        u_max_i = std::max(u_max_i, u_view_o.get(c));
                        u_min_i = std::min(u_min_i, u_view_o.get(c));
                    }
                }

                Q_plus_view.set(r, (Mii / dt) * (u_max_i - u_view_d.get(r)));
                Q_minus_view.set(r, (Mii / dt) * (u_min_i - u_view_d.get(r)));
            }
        }

        // Eq (66) + (82)
        static void create_update(Matrix &F,
                                  const Vector &R_minus,
                                  const Vector &R_minus_ghosts,
                                  const Vector &R_plus,
                                  const Vector &R_plus_ghosts,
                                  Vector &f_bar) {
            auto R_plus_view_d = local_view_host(R_plus);
            auto R_minus_view_d = local_view_host(R_minus);

            auto R_plus_view_o = local_view_host(R_plus_ghosts);
            auto R_minus_view_o = local_view_host(R_minus_ghosts);

            auto f_bar_view = local_view_host(f_bar);

            PetscCrsView F_d, F_o;
            views_host(F, F_d, F_o);

            SizeType rows = F_d.rows();
            for (SizeType r = 0; r < rows; r++) {
                Scalar fi = 0;

                {
                    auto row = F_d.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        const SizeType c = row.colidx(k);
                        if (c == r) continue;

                        Scalar Fij = row.value(k);
                        Scalar alpha_ij = 0;

                        if (Fij > 0) {
                            alpha_ij = std::min(R_plus_view_d.get(r), R_minus_view_d.get(c));
                        } else {
                            alpha_ij = std::min(R_minus_view_d.get(r), R_plus_view_d.get(c));
                        }

                        fi += alpha_ij * Fij;
                    }
                }

                if(!empty(F_o))
                {
                    auto row = F_o.row(r);
                    for (SizeType k = 0; k < row.length; k++) {
                        const SizeType c = row.colidx(k);

                        Scalar Fij = row.value(k);

                        Scalar alpha_ij = 0;

                        if (Fij > 0) {
                            alpha_ij = std::min(R_plus_view_d.get(r), R_minus_view_o.get(c));
                        } else {
                            alpha_ij = std::min(R_minus_view_d.get(r), R_plus_view_o.get(c));
                        }

                        fi += alpha_ij * Fij;
                    }
                }

                f_bar_view.set(r, fi);
            }
        }

        void post_process(const Vector &g, Vector &x) const override {
            UTOPIA_TRACE_SCOPE("StabilizeTransport::post_process");
            //!!! -g instead of g because the rhs is = -g
            Vector x_dot = (*A_corrected_) * x - g;
            // Vector x_dot = (*A_corrected_) * x + g;
            x_dot = e_mul(inverse_mass_vector_, x_dot);

            Matrix F = *A_corrected_;
            F *= 0;

            IndexArray ghosts;
            A_corrected_->ghosts(ghosts);

            // Split vectors into local block and off diagonal block
            // retrieve ghost values of vectors
            Vector x_ghosts, x_dot_ghosts;
            x.select(ghosts, x_ghosts);
            x_dot.select(ghosts, x_dot_ghosts);

            add_contrib(*A_diff_, x, x_ghosts, F);
            add_contrib(*M_diff_, x, x_ghosts, F);

            Vector P_plus(layout(x), 0), P_minus(layout(x), 0);
            negative_postive_antidiffusive_fluxes(F, P_plus, P_minus);

            Vector Q_plus(layout(x), 0), Q_minus(layout(x), 0);
            min_max_bound(dt_, *M_corrected_, x, x_ghosts, Q_minus, Q_plus);

            // Book_Kuzmin.pdf page 163 eq (82)
            Vector R_plus(layout(x), 0), R_minus(layout(x), 0);
            R_plus = Q_plus / P_plus;
            R_plus.e_min(1);

            R_minus = Q_minus / P_minus;
            R_minus.e_min(1);

            Vector R_plus_ghosts, R_minus_ghosts;
            R_plus.select(ghosts, R_plus_ghosts);
            R_minus.select(ghosts, R_minus_ghosts);

            Vector f_bar(layout(x), 0);
            create_update(F, R_minus, R_minus_ghosts, R_plus, R_minus_ghosts, f_bar);

            f_bar *= dt_;
            f_bar = e_mul(inverse_mass_vector_, f_bar);
            x += f_bar;
        }

    private:
        std::shared_ptr<Matrix> A_corrected_;
        std::shared_ptr<Matrix> A_diff_;
        std::shared_ptr<Matrix> M_corrected_;
        std::shared_ptr<Matrix> M_diff_;
        Vector inverse_mass_vector_;
        Scalar dt_{1};
    };

    template <class Matrix>
    void register_transfomers() {
        static bool once = false;

        if (!once) {
            MatrixTransformer<Matrix>::template Register<StabilizeTransport<Matrix>>("StabilizeTransport");
            once = true;
        }
    }

}  // namespace utopia

#endif  // UTOPIA_MATRIX_TRANSFORMER_HPP

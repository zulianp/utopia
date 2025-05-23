#ifndef UTOPIA_LIBMESH_TRANSFEROPERATOR_HPP
#define UTOPIA_LIBMESH_TRANSFEROPERATOR_HPP

#include "utopia.hpp"

#include "utopia_Path.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_Local2Global.hpp"
#include "utopia_LocalAssembler.hpp"

#include <memory>

namespace utopia {

    class TransferOptions : public Configurable {
    public:
        TransferOptions()
            : from_var_num(0), to_var_num(0), n_var(1), to_trace_space(false), tags({}), has_covering(true) {}

        void read(Input &in) override {
            in.get("from-var-num", from_var_num);
            in.get("to-var-num", to_var_num);
            in.get("n-var", n_var);
            in.get("has_covering", has_covering);

            in.get("coupling", [this](Input &in) {
                in.get_all([this](Input &in) {
                    int master = -1, slave = -1;

                    in.get("master", master);
                    in.get("slave", slave);

                    assert(master != -1);
                    assert(slave != -1);

                    tags.emplace_back(master, slave);
                });
            });
        }

        int from_var_num;
        int to_var_num;
        int n_var;
        bool to_trace_space;
        std::vector<std::pair<int, int>> tags;
        bool has_covering;
    };

    class TransferOperator {
    public:
        virtual ~TransferOperator() {}
        virtual void apply(const UVector &from, UVector &to) const = 0;
        virtual void apply_transpose(const UVector &from, UVector &to) const = 0;
        virtual void describe(std::ostream &) const {}
        virtual bool write(const Path &) const { return false; }
    };

    /**
     * @brief constructed as (D^-1 * B) * ( . )
     */
    class L2TransferOperator final : public TransferOperator {
    public:
        inline void apply(const UVector &from, UVector &to) const override {
            UVector B_from = *B * from;

            if (empty(to)) {
                to = B_from;
            }

            linear_solver->apply(B_from, to);
        }

        ///@brief assumes that D is symmetric
        void apply_transpose(const UVector &from, UVector &to) const override {
            UVector D_inv_from = local_zeros(local_size(*D).get(0));
            linear_solver->apply(from, D_inv_from);
            to = transpose(*B) * D_inv_from;
        }

        static std::unique_ptr<L2TransferOperator> Make(
            const std::shared_ptr<USparseMatrix> &B,
            const std::shared_ptr<USparseMatrix> &D,
            const std::shared_ptr<LinearSolver<USparseMatrix, UVector>> &linear_solver =
                std::make_shared<BiCGStab<USparseMatrix, UVector>>(),
            const double tol = 1e-16) {
            auto ret = utopia::make_unique<L2TransferOperator>(B, D, linear_solver);
            ret->init();
            return ret;
        }

        static std::unique_ptr<L2TransferOperator> make_fixed(
            const std::shared_ptr<USparseMatrix> &B,
            const std::shared_ptr<USparseMatrix> &D,
            const std::shared_ptr<LinearSolver<USparseMatrix, UVector>> &linear_solver =
                std::make_shared<BiCGStab<USparseMatrix, UVector>>(),
            const double tol = 1e-16) {
            auto ret = utopia::make_unique<L2TransferOperator>(B, D, linear_solver);
            ret->fix_mass_matrix_operator(tol);
            ret->init();
            return ret;
        }

        static std::unique_ptr<L2TransferOperator> make_restricted(
            const std::shared_ptr<USparseMatrix> &B,
            const std::shared_ptr<USparseMatrix> &D,
            const std::shared_ptr<LinearSolver<USparseMatrix, UVector>> &linear_solver =
                std::make_shared<BiCGStab<USparseMatrix, UVector>>(),
            const double tol = 1e-16) {
            auto ret = utopia::make_unique<L2TransferOperator>(B, D, linear_solver);
            ret->restrict_mass_matrix(tol);
            ret->init();
            return ret;
        }

        inline void describe(std::ostream &os) const override {
            UVector t_from = local_values(local_size(*B).get(1), 1);
            UVector t_to;
            apply(t_from, t_to);

            double t_max = max(t_to);
            double t_min = min(t_to);

            double sum_D = sum(*D);
            double sum_B = sum(*B);

            os << "------------------------------------------\n";
            os << "L2TransferOperator:\n";
            os << "row sum [" << t_min << ", " << t_max << "] subset of [0, 1]" << std::endl;
            os << "sum(B) = " << sum_B << ", sum(D) = " << sum_D << std::endl;
            os << "------------------------------------------\n";
        }

        bool write(const Path &path) const override {
            return utopia::write(path / "B.m", *B) && utopia::write(path / "D.m", *D);
        }

        inline L2TransferOperator(const std::shared_ptr<USparseMatrix> &B,
                                  const std::shared_ptr<USparseMatrix> &D,
                                  const std::shared_ptr<LinearSolver<USparseMatrix, UVector>> &linear_solver =
                                      std::make_shared<BiCGStab<USparseMatrix, UVector>>())
            : B(B), D(D), linear_solver(linear_solver) {
            assert(B);
            assert(D);
            assert(linear_solver);
        }

        inline void init() { linear_solver->update(D); }

        void restrict_mass_matrix_old(const double tol = 1e-16) {
            auto rr = row_range(*B);

            const SizeType n_local = rr.extent();

            std::vector<bool> exists(n_local, false);
            std::vector<USparseMatrix::SizeType> indices;
            indices.reserve(n_local);

            utopia::each_read(*B, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
                if (std::abs(value) > tol) {
                    auto idx = i - rr.begin();
                    exists[idx] = true;
                }
            });

            for (SizeType i = 0; i < n_local; ++i) {
                if (!exists[i]) {
                    indices.push_back(rr.begin() + i);
                }
            }

            set_zero_rows(*D, indices, 1.);
        }

        void restrict_mass_matrix(const double tol = 1e-16) {
            auto rr = row_range(*B);

            const SizeType n_local = rr.extent();

            std::vector<bool> exists(n_local, false);
            std::vector<USparseMatrix::SizeType> indices;
            indices.reserve(n_local);

            UVector sum_D = sum(*D, 1);
            UVector sum_B = sum(*B, 1);

            Read<UVector> r_D(sum_D), r_B(sum_B);

            for (auto i = rr.begin(); i < rr.end(); ++i) {
                if (!approxeq(sum_D.get(i), sum_B.get(i), tol)) {
                    indices.push_back(rr.begin() + i);
                }
            }

            set_zero_rows(*D, indices, 1.);
        }

        void fix_mass_matrix_operator(const double tol = 1e-16) {
            UVector d;

            Size s = local_size(*D);
            d = local_values(s.get(0), 1.);

            {
                Write<UVector> w_d(d);

                each_read(*D, [&d, tol](const SizeType i, const SizeType, const double val) {
                    if (std::abs(val) > tol) {
                        d.set(i, 0.);
                    }
                });
            }

            (*D) += USparseMatrix(diag(d));
        }

        friend Size local_size(const L2TransferOperator &op) { return local_size(*op.B); }

    private:
        std::shared_ptr<USparseMatrix> B;
        std::shared_ptr<USparseMatrix> D;
        std::shared_ptr<LinearSolver<USparseMatrix, UVector>> linear_solver;
    };

    class PseudoL2TransferOperator final : public TransferOperator {
    public:
        inline void apply(const UVector &from, UVector &to) const override {
            assert(T);
            to = *T * from;
        }

        inline void apply_transpose(const UVector &from, UVector &to) const override {
            assert(T);
            to = transpose(*T) * from;
        }

        PseudoL2TransferOperator() {}

        inline void init_from_coupling_operator(const USparseMatrix &B) {
            T = std::make_shared<USparseMatrix>();
            UVector d = sum(B, 1);

            {
                ReadAndWrite<UVector> rw_(d);
                auto r = range(d);
                for (auto k = r.begin(); k != r.end(); ++k) {
                    if (approxeq(d.get(k), 0.0, 1e-14)) {
                        d.set(k, 1.);
                    }
                }
            }

            *T = diag(1. / d) * B;
        }

        inline void init_from_coupling_and_mass_operator(const USparseMatrix &B, const USparseMatrix &M) {
            T = std::make_shared<USparseMatrix>();
            UVector d = sum(M, 1);

            {
                ReadAndWrite<UVector> rw_(d);
                auto r = range(d);
                for (auto k = r.begin(); k != r.end(); ++k) {
                    if (approxeq(d.get(k), 0.0, 1e-14)) {
                        d.set(k, 1.);
                    }
                }
            }

            *T = diag(1. / d) * B;
        }

        PseudoL2TransferOperator(const std::shared_ptr<USparseMatrix> &T) : T(T) { assert(T); }

        inline void describe(std::ostream &os) const override {
            UVector t = sum(*T, 1);
            double t_max = max(t);
            double t_min = min(t);
            double t_sum = sum(t);

            os << "------------------------------------------\n";
            os << "PseudoL2TransferOperator:\n";
            os << "row sum [" << t_min << ", " << t_max << "] subset of [0, 1]" << std::endl;
            os << "sum(T): " << t_sum << " <= " << size(*T).get(0) << "\n";
            os << "------------------------------------------\n";
        }

        bool write(const Path &path) const override { return utopia::write(path / "T.m", *T); }

        inline std::shared_ptr<USparseMatrix> matrix() { return T; }

    private:
        std::shared_ptr<USparseMatrix> T;
    };

    class PermutedOperator final : public TransferOperator {
    public:
        PermutedOperator(const std::shared_ptr<TransferOperator> &op,
                         const std::shared_ptr<USparseMatrix> &from_permutation,
                         const std::shared_ptr<USparseMatrix> &to_permutation)
            : op_(op), from_permutation_(from_permutation), to_permutation_(to_permutation) {
            from_buffer_ = utopia::make_unique<UVector>();
            to_buffer_ = utopia::make_unique<UVector>();
        }

        inline void apply(const UVector &from, UVector &to) const {
            if (from_permutation_) {
                *from_buffer_ = *from_permutation_ * from;
            } else {
                *from_buffer_ = from;
            }

            op_->apply(*from_buffer_, *to_buffer_);

            if (to_permutation_) {
                to = transpose(*to_permutation_) * (*to_buffer_);
            } else {
                to = *to_buffer_;
            }
        }

        inline void apply_transpose(const UVector &to, UVector &from) const {
            if (to_permutation_) {
                *to_buffer_ = *to_permutation_ * to;
            } else {
                *to_buffer_ = to;
            }

            op_->apply_transpose(*to_buffer_, *from_buffer_);

            if (from_permutation_) {
                from = transpose(*from_permutation_) * (*from_buffer_);
            } else {
                from = *from_buffer_;
            }
        }

        std::shared_ptr<TransferOperator> op_;
        std::shared_ptr<USparseMatrix> from_permutation_;
        std::shared_ptr<USparseMatrix> to_permutation_;

        std::unique_ptr<UVector> from_buffer_, to_buffer_;
    };

    class NormalizedOperator final : public TransferOperator {
    public:
        using Scalar = UTOPIA_SCALAR(UVector);

        NormalizedOperator(const Size &op_local_size, const std::shared_ptr<TransferOperator> &op) : op_(op) {
            UVector ones = local_values(op_local_size.get(1), 1.);
            op->apply(ones, rescale_);
            inv_rescale_ = 1. / rescale_;
            temp_ = utopia::make_unique<UVector>();
        }

        inline void apply(const UVector &from, UVector &to) const {
            op_->apply(from, *temp_);
            to = e_mul(inv_rescale_, *temp_);
        }

        inline void apply_transpose(const UVector &from, UVector &to) const {
            *temp_ = e_mul(rescale_, from);
            op_->apply_transpose(*temp_, to);
        }

        inline void describe(std::ostream &os) const {
            os << "non normalized: \n";
            op_->describe(os);
        }

    private:
        std::shared_ptr<TransferOperator> op_;
        UVector rescale_, inv_rescale_;
        std::unique_ptr<UVector> temp_;
    };

    class ClampedOperator final : public TransferOperator {
    public:
        using Scalar = UTOPIA_SCALAR(UVector);

        ClampedOperator(const std::shared_ptr<TransferOperator> &op) : op_(op) {}

        inline void apply(const UVector &from, UVector &to) const {
            op_->apply(from, to);
            clamp(from, to);
        }

        inline void apply_transpose(const UVector &from, UVector &to) const {
            op_->apply_transpose(from, to);
            clamp(from, to);
        }

    private:
        static void clamp(const UVector &from, UVector &to) {
            const Scalar min_val = min(from);
            const Scalar max_val = max(from);

            ReadAndWrite<UVector> rw_(to);

            auto r = range(to);
            for (auto i = r.begin(); i < r.end(); ++i) {
                const Scalar val = to.get(i);
                to.set(i, std::min(max_val, std::max(min_val, val)));
            }
        }

        std::shared_ptr<TransferOperator> op_;
    };

    class ForceZeroExtension final : public TransferOperator {
    public:
        using Scalar = UTOPIA_SCALAR(UVector);

        ForceZeroExtension(const std::shared_ptr<TransferOperator> &op, const double tol) : op_(op), tol_(tol) {}

        inline void apply(const UVector &from, UVector &to) const {
            op_->apply(from, to);
            clamp(from, to);
        }

        inline void apply_transpose(const UVector &from, UVector &to) const {
            op_->apply_transpose(from, to);
            clamp(from, to);
        }

    private:
        void clamp(const UVector &from, UVector &to) const {
            const Scalar min_val = min(from);
            const Scalar max_val = max(from);

            ReadAndWrite<UVector> rw_(to);

            auto r = range(to);
            for (auto i = r.begin(); i < r.end(); ++i) {
                Scalar val = to.get(i);

                if (val + tol_ < min_val) {
                    to.set(i, 0.);
                } else if (val - tol_ > max_val) {
                    to.set(i, 0.);
                } else {
                    to.set(i, val);
                }
            }
        }

        std::shared_ptr<TransferOperator> op_;
        double tol_;
    };

    class BidirectionalOperator final : public TransferOperator {
    public:
        BidirectionalOperator(const std::shared_ptr<TransferOperator> &forward,
                              const std::shared_ptr<TransferOperator> &backward)
            : forward_(forward), backward_(backward) {}

        inline void apply(const UVector &from, UVector &to) const { forward_->apply(from, to); }

        inline void apply_transpose(const UVector &from, UVector &to) const { backward_->apply(from, to); }

        inline void describe(std::ostream &os) const {
            if (forward_) forward_->describe(os);
            if (backward_) backward_->describe(os);
        }

        inline bool write(const Path &) const { return false; }

    private:
        std::shared_ptr<TransferOperator> forward_;
        std::shared_ptr<TransferOperator> backward_;
    };

    class Interpolator final : public TransferOperator {
    public:
        inline void apply(const UVector &from, UVector &to) const override { to = *T * from; }

        void apply_transpose(const UVector &from, UVector &to) const override {
            assert(T);
            to = transpose(*T) * from;
        }

        Interpolator(const std::shared_ptr<USparseMatrix> &T) : T(T) { assert(T); }

        void normalize_rows() {
            UVector d = sum(*T, 1);

            auto r = range(d);

            {
                ReadAndWrite<UVector> rw_(d);
                for (auto k = r.begin(); k != r.end(); ++k) {
                    if (approxeq(d.get(k), 0.0, 1e-14)) {
                        d.set(k, 1.);
                    }
                }
            }

            *T = diag(1. / d) * (*T);
        }

        inline void describe(std::ostream &os) const override {
            UVector t = sum(*T, 1);
            double t_max = max(t);
            double t_min = min(t);
            double t_sum = sum(t);

            os << "------------------------------------------\n";
            os << "Interpolator:\n";
            os << "row sum [" << t_min << ", " << t_max << "] subset of [0, 1]" << std::endl;
            os << "sum(T): " << t_sum << " <= " << size(*T).get(0) << "\n";
            os << "------------------------------------------\n";
        }

        bool write(const Path &path) const override { return utopia::write(path / "T.m", *T); }

        inline std::shared_ptr<USparseMatrix> matrix() { return T; }

    private:
        std::shared_ptr<USparseMatrix> T;
    };

    enum TransferOperatorType {
        INTERPOLATION = 0,
        L2_PROJECTION = 1,
        PSEUDO_L2_PROJECTION = 2,
        APPROX_L2_PROJECTION = 3,
        BIDIRECTIONAL_L2_PROJECTION = 4,
        BIDIRECTIONAL_PSEUDO_L2_PROJECTION = 5,
    };
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_TRANSFEROPERATOR_HPP

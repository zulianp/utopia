#ifndef UTOPIA_KOKKOS_OP_HPP
#define UTOPIA_KOKKOS_OP_HPP

#include "utopia_Base.hpp"

namespace utopia {

    class TestOp {
    public:
        UTOPIA_INLINE_FUNCTION constexpr TestOp(const int offset_test = 0) : offset_test_(offset_test) {}
        UTOPIA_INLINE_FUNCTION constexpr int offset_test() const { return offset_test_; }
        UTOPIA_INLINE_FUNCTION void set_offset_test(const int offset_test) { offset_test_ = offset_test; }
        int offset_test_{0};
    };

    class TrialOp {
    public:
        UTOPIA_INLINE_FUNCTION constexpr TrialOp(const int offset_trial = 0) : offset_trial_(offset_trial) {}
        UTOPIA_INLINE_FUNCTION constexpr int offset_trial() const { return offset_trial_; }
        UTOPIA_INLINE_FUNCTION void set_offset_trial(const int offset_trial) { offset_trial_ = offset_trial; }
        int offset_trial_{0};
    };

    class TestTrialOp : public TestOp, public TrialOp {
    public:
        UTOPIA_INLINE_FUNCTION constexpr TestTrialOp(const int offset_test = 0, const int offset_trial = 0)
            : TestOp(offset_test), TrialOp(offset_trial) {}
    };

    class TestOpParams : public TestOp, public Configurable {
    public:
        void read(Input &in) override { in.get("offset_test", offset_test_); }

        template <class Op>
        inline Op update(Op &&op) const {
            op.set_offset_test(offset_test());
            return std::forward<Op>(op);
        }
    };

    class TrialOpParams : public TrialOp, public Configurable {
    public:
        void read(Input &in) override { in.get("offset_trial", offset_trial_); }

        template <class Op>
        inline Op update(Op &&op) const {
            op.set_offset_trial(offset_trial());
            return std::forward<Op>(op);
        }
    };

    class TestTrialOpParams : public TestTrialOp, public Configurable {
    public:
        void read(Input &in) override {
            in.get("offset_test", offset_test_);
            in.get("offset_trial", offset_trial_);
        }

        template <class Op>
        inline Op update(Op &&op) const {
            op.set_offset_test(offset_test());
            op.set_offset_trial(offset_trial());
            return std::forward<Op>(op);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_KOKKOS_OP_HPP

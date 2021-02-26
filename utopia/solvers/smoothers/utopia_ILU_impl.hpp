#ifndef UTOPIA_ILU_IMPL_HPP
#define UTOPIA_ILU_IMPL_HPP

#include "utopia_DILUDecompose.hpp"
#include "utopia_ILU.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Backend>
    class ILU<Matrix, Vector, Backend>::Impl final : public Configurable {
    public:
        // using DefaultILU_t = utopia::DILUAlgorithm<Matrix, Vector>;
        using DefaultILU_t = utopia::ILUDecompose<Matrix>;

        void update(const Matrix &mat) { algo_->update(mat); }

        void apply(const Vector &in, Vector &out) { algo_->apply(in, out); }

        void read(Input &in) override {
            int block_size = 0;
            in.get("block_size", block_size);

            bool use_dilu = false;

            in.get("use_dilu", use_dilu);

            if (use_dilu) {
                algo_ = utopia::make_unique<DILUAlgorithm<Matrix, Vector>>();
            } else {
                // if (block_size == 1) {
                //     algo_ = utopia::make_unique<BlockILUAlgorithm<Matrix, 1>>();
                // } else
                if (block_size == 2) {
                    algo_ = utopia::make_unique<BlockILUAlgorithm<Matrix, 2>>();
                } else if (block_size == 3) {
                    algo_ = utopia::make_unique<BlockILUAlgorithm<Matrix, 3>>();
                } else if (block_size == 4) {
                    algo_ = utopia::make_unique<BlockILUAlgorithm<Matrix, 4>>();
                }
                //  else {
                //     algo_ = utopia::make_unique<DefaultILU_t>();
                // }
            }

            algo_->read(in);
        }

        Impl() : algo_(utopia::make_unique<DefaultILU_t>()) {}

        Matrix decomposition;
        Vector residual, correction;
        std::unique_ptr<ILUAlgorithm<Matrix, Vector>> algo_;
    };

    template <class Matrix, class Vector, int Backend>
    ILU<Matrix, Vector, Backend>::~ILU() = default;

    template <class Matrix, class Vector, int Backend>
    ILU<Matrix, Vector, Backend>::ILU() : impl_(utopia::make_unique<Impl>()) {}

    template <class Matrix, class Vector, int Backend>
    ILU<Matrix, Vector, Backend> *ILU<Matrix, Vector, Backend>::clone() const {
        return new ILU();
    }

    template <class Matrix, class Vector, int Backend>
    void ILU<Matrix, Vector, Backend>::read(Input &in) {
        Super::read(in);
        impl_->read(in);
    }

    template <class Matrix, class Vector, int Backend>
    void ILU<Matrix, Vector, Backend>::print_usage(std::ostream &os) const {
        Super::print_usage(os);
    }

    template <class Matrix, class Vector, int Backend>
    bool ILU<Matrix, Vector, Backend>::smooth(const Vector &b, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("ILU::smooth");
        auto &A = *this->get_operator();

        SizeType n_sweeps = this->sweeps();

        for (int i = 0; i < n_sweeps; ++i) {
            A.apply(x, impl_->residual);
            impl_->residual = b - impl_->residual;
            impl_->apply(impl_->residual, impl_->correction);
            x += impl_->correction;
        }

        UTOPIA_TRACE_REGION_END("ILU::smooth");
        return true;
    }

    template <class Matrix, class Vector, int Backend>
    bool ILU<Matrix, Vector, Backend>::apply(const Vector &b, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("ILU::apply");

        if (this->verbose()) {
            this->init_solver("ILU", {" it. ", "|| residual ||"});
        }

        auto &A = *this->get_operator();

        bool converged = false;

        int iteration = 0;
        while (!converged) {
            // smooth(b, x);

            A.apply(x, impl_->residual);
            impl_->residual = b - impl_->residual;
            impl_->apply(impl_->residual, impl_->correction);
            x += impl_->correction;

            const Scalar norm_r = norm2(impl_->residual);

            if (this->verbose()) {
                PrintInfo::print_iter_status(iteration, {norm_r});
            }

            converged = this->check_convergence(iteration, 1, 1, norm_r);
            ++iteration;

            if (converged) break;
        }

        UTOPIA_TRACE_REGION_END("ILU::apply");
        return converged;
    }

    template <class Matrix, class Vector, int Backend>
    void ILU<Matrix, Vector, Backend>::init_memory(const Layout &layout) {
        if (impl_->residual.empty() || !layout.same(utopia::layout(impl_->residual))) {
            impl_->residual.zeros(layout);
            impl_->correction.zeros(layout);
        }
    }

    template <class Matrix, class Vector, int Backend>
    void ILU<Matrix, Vector, Backend>::update(const std::shared_ptr<const Matrix> &op) {
        UTOPIA_TRACE_REGION_BEGIN("ILU::update");

        Super::update(op);
        impl_->update(*op);
        init_memory(row_layout(*op));

        UTOPIA_TRACE_REGION_END("ILU::update");
    }

}  // namespace utopia

#endif  // UTOPIA_ILU_IMPL_HPP

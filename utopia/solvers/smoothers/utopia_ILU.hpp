#ifndef UTOPIA_ILU_HPP
#define UTOPIA_ILU_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

#include "utopia_IterativeSolver.hpp"

#include <memory>

namespace utopia {

    template <class Matrix, class Vector>
    class ILUAlgorithm : public Configurable, public Clonable {
    public:
        virtual ~ILUAlgorithm() = default;
        virtual bool update(const Matrix &mat) = 0;
        virtual void apply(const Vector &b, Vector &x) = 0;
        void read(Input &) override {}
        ILUAlgorithm *clone() const override = 0;
    };

    template <class Matrix, int BlockSize>
    class BlockILUAlgorithm {};

    template <class Matrix, int Backend = Traits<Matrix>::Backend>
    class ILUDecompose {};

    template <class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class ILU final : public IterativeSolver<Matrix, Vector> {
    public:
        using Layout = typename Traits<Vector>::Layout;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Super = utopia::IterativeSolver<Matrix, Vector>;

        ~ILU() override;
        ILU();

        ILU *clone() const override;
        void read(Input &in) override;
        void print_usage(std::ostream &os) const override;
        bool smooth(const Vector &b, Vector &x) override;
        bool apply(const Vector &b, Vector &x) override;
        void init_memory(const Layout &layout) override;
        void update(const std::shared_ptr<const Matrix> &op) override;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_ILU_HPP

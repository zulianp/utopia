#ifndef UTOPIA_TPETRA_OPERATOR_HPP
#define UTOPIA_TPETRA_OPERATOR_HPP

#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_Operator.hpp>

#include "utopia_Tpetra_Matrix.hpp"

#include <cassert>

#include "utopia_Instance.hpp"

namespace utopia {

    template <class Vector>
    class TpetraOperator : public Tpetra::Operator<typename Traits<Vector>::Scalar,
                                                   typename Traits<Vector>::LocalSizeType,
                                                   typename Traits<Vector>::SizeType,
                                                   typename Traits<Vector>::Node> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using LocalSizeType = typename Traits::LocalSizeType;
        using SizeType = typename Traits::SizeType;
        using Node = typename Traits::Node;
        using MapType = Tpetra::Map<LocalSizeType, SizeType, Node>;
        using MultiVectorType = Tpetra::MultiVector<Scalar, LocalSizeType, SizeType, Node>;
        using VectorType = Tpetra::Vector<Scalar, LocalSizeType, SizeType, Node>;

        TpetraOperator(const std::shared_ptr<Operator<Vector>> &op) : op_(op) { assert(op); }

        static void wrap_mv(MultiVectorType &in, Vector &out) {
            Teuchos::RCP<VectorType> yy;
            auto sub_ptr = dynamic_cast<VectorType *>(&in);
            if (!sub_ptr) {
                assert(in.getNumVectors() == 1);
                yy = Teuchos::RCP<VectorType>(new VectorType(in, 0), true);
            } else {
                yy = Teuchos::RCP<VectorType>(sub_ptr, false);
            }

            out.wrap(yy);
        }

        void apply(const MultiVectorType &X,
                   MultiVectorType &Y,
                   Teuchos::ETransp mode = Teuchos::NO_TRANS,
                   Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                   Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const override {
            Vector wx, wy;

            wrap_mv(const_cast<MultiVectorType &>(X), wx);
            wrap_mv(Y, wy);

            op_->apply(wx, wy);

            wx.unwrap(Teuchos::RCP<VectorType>());
            wy.unwrap(Teuchos::RCP<VectorType>());
        }

        Teuchos::RCP<const Tpetra::Map<LocalSizeType, SizeType, Node>> getDomainMap() const override {
            auto lsz = op_->local_size();
            auto gsz = op_->size();
            auto comm = op_->comm().get();
            return Teuchos::rcp(new MapType(gsz.get(1), lsz.get(1), 0, comm));
        }

        Teuchos::RCP<const Tpetra::Map<LocalSizeType, SizeType, Node>> getRangeMap() const override {
            auto lsz = op_->local_size();
            auto gsz = op_->size();
            auto comm = op_->comm().get();
            return Teuchos::rcp(new MapType(gsz.get(0), lsz.get(0), 0, comm));
        }

        bool hasTransposeApply() const override { return false; }
        bool hasDiagonal() const override { return false; }
        void getLocalDiagCopy(Tpetra::Vector<Scalar, LocalSizeType, SizeType, Node> & /*diag*/) const {
            assert(false);
            Utopia::Abort("TpetraOperator: IMPLEMENT ME!");
        }

    private:
        std::shared_ptr<Operator<Vector>> op_;
    };

    template <class Vector>
    inline Teuchos::RCP<TpetraOperator<Vector>> tpetra_operator(const std::shared_ptr<Operator<Vector>> &op) {
        return Teuchos::rcp(new TpetraOperator<Vector>(op));
    }

}  // namespace utopia

#endif  // UTOPIA_TPETRA_OPERATOR_HPP

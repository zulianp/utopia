#ifndef UTOPIA_CONSTRAINT_QP_PMPRGP
#define UTOPIA_CONSTRAINT_QP_PMPRGP

#include <string>
#include "utopia_Algorithms.hpp"
#include "utopia_Allocations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_For.hpp"
#include "utopia_QPSolver.hpp"

#include "utopia_MPRGP.hpp"

#include "utopia_petsc_ProjectedGaussSeidel.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class PMPRGP final : public MPRGP<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using Solver = utopia::LinearSolver<Matrix, Vector>;
        using Super = utopia::MPRGP<Matrix, Vector>;

    public:
        using Super::solve;
        using Super::update;

        PMPRGP() {  MPRGP<Matrix, Vector>::set_preconditioner(pgs()); }

        static std::shared_ptr<ProjectedGaussSeidel<Matrix, Vector>> pgs() {
            auto ret = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            ret->max_it(1);
            ret->n_local_sweeps(4);
            return ret;
        }

        void read(Input &in) override {
            MPRGP<Matrix, Vector>::read(in);


            std::string preconditioner_type = "pgs";
            in.require("preconditioner_type", preconditioner_type);

            std::shared_ptr<Preconditioner<Vector>> prec;
            if (preconditioner_type == "pgs") {
                prec = pgs();
            } else {
                Utopia::Abort("Only pgs is supported right now!");
            }

            in.get("preconditioner", *prec);
            MPRGP<Matrix, Vector>::set_preconditioner(prec);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_CONSTRAINT_QP_PMPRGP

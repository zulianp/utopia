#ifndef UTOPIA_M3ELINSOL_HPP
#define UTOPIA_M3ELINSOL_HPP

#include "utopia_IterativeSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_Input.hpp"

namespace utopia {
    class M3ELinSol {
    public:
        class Impl;
        ~M3ELinSol();
        M3ELinSol();

        std::shared_ptr<Impl> impl;
    };

    template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class ASPAMG final : public IterativeSolver<Matrix, Vector> {
    public:
        bool apply(const Vector &b, Vector &x) override;

        /*! @brief if overriden the subclass has to also call this one first
         */
        virtual void update(const std::shared_ptr<const Matrix> &op) override;

        void print_system(const bool binwrite, const std::string systemfile);

        ASPAMG * clone() const override
        {
            //FIXME
            return new ASPAMG();
        }

        void read(Input &in) override;


        ASPAMG()
        {}

    private:
        M3ELinSol solver;

    };
}

#endif //UTOPIA_M3ELINSOL_HPP

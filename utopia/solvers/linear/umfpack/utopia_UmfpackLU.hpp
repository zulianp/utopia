
#ifndef UTOPIA_UTOPIA_UMFPACKLU_HPP
#define UTOPIA_UTOPIA_UMFPACKLU_HPP

#include "utopia_DirectSolver.hpp"
#include "utopia_BLAS.hpp"
#include "utopia_LinearSolverInterfaces.hpp"

namespace utopia {

  namespace internals {
    class UmfpackLU : virtual public DirectSolver<CRSMatrixd, Vectord>,
                      virtual public DirectSolver<CCSMatrixd, Vectord> {
    public:

        UmfpackLU()
        {}

        ~UmfpackLU();

        bool solve(const CCSMatrixd &mat, const Vectord &rhs, Vectord &solution);
        bool solve(const CRSMatrixd &matrix, const Vectord &rhs, Vectord &solution);

        inline const double * ptr(const Vectord &v) const
        {
            return &v.implementation()[0];
        }

        inline double * ptr(Vectord &v) const
        {
            return &v.implementation()[0];
        }

        inline const double * ptr(const CRSMatrixd &m) const
        {
            return &m.implementation().getEntries()[0];
        }

        virtual bool apply(const Vectord &rhs, Vectord &sol)
        {
           //deal with multiple inheritance
           if(DirectSolver<CRSMatrixd, Vectord>::has_operator()) {
               return DirectSolver<CRSMatrixd, Vectord>::apply(rhs, sol);
           } else {
               return DirectSolver<CCSMatrixd, Vectord>::apply(rhs, sol);
           }
        }

    private:

        class Internal {
        public:
          void * symbolic, * numeric;
          std::vector<double> control;
          std::vector<double> info;

            Internal()
                    : symbolic(NULL), numeric(NULL)
            {}
        };

        bool init(const CCSMatrixd &mat, Internal &internal) const;
        void cleanUp(Internal &internal) const;

        bool solve(const CCSMatrixd &mat, Internal &internal, const double * rhs, double * solution);
    };
  }

    template<>
    class LUDecomposition<CRSMatrixd, Vectord, BLAS> : public DirectSolver<CRSMatrixd, Vectord> {
    public:
      inline bool apply(const Vectord &rhs, Vectord &sol)
      {
         return strategy_.apply(rhs, sol);
      }

    private:
      internals::UmfpackLU strategy_;
    };

    template<>
    class LUDecomposition<CCSMatrixd, Vectord, BLAS> : public DirectSolver<CCSMatrixd, Vectord> {
    public:
      inline bool apply(const Vectord &rhs, Vectord &sol)
      {
         return strategy_.apply(rhs, sol);
      }

    private:
      internals::UmfpackLU strategy_;
    };
}

#endif //UTOPIA_UTOPIA_UMFPACKLU_HPP

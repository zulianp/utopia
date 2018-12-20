#ifndef UTOPIA_BICG_STAB_HPP
#define UTOPIA_BICG_STAB_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Size.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_PreconditionedSolver.hpp"

namespace utopia {

	template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
	class BiCGStab : public PreconditionedSolver<Matrix, Vector>, public Smoother<Matrix, Vector>, public MatrixFreeLinearSolver<Vector> {
	public:
		typedef UTOPIA_SCALAR(Vector) 	 Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
		typedef utopia::Preconditioner<Vector> Preconditioner;

		using PreconditionedSolver<Matrix, Vector>::solve;
	
		BiCGStab();	
		BiCGStab * clone() const override;
		
		inline bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override
		{
			if(this->get_preconditioner()) {
				return solve_preconditioned(A, b, x);
			} else {
				return solve_unpreconditioned(A, b, x);
			}
		}

		inline bool apply(const Vector &b, Vector &x) override
		{
			auto A_ptr = utopia::op(this->get_operator());
			return solve(*A_ptr, b, x);
		}

		bool smooth(const Vector &rhs, Vector &x) override;

		//for chosing the preconditioned solver one
		void update(const std::shared_ptr<const Matrix> &op) override;


		void read(Input &in) override
        {
            MatrixFreeLinearSolver<Vector>::read(in);
            // Smoother<Matrix, Vector>::read(in); 
            // PreconditionedSolver<Matrix, Vector>::read(in); 
        }


        void print_usage(std::ostream &os) const override
        {
            MatrixFreeLinearSolver<Vector>::print_usage(os);
            // Smoother<Matrix, Vector>::print_usage(os);
            // PreconditionedSolver<Matrix, Vector>::print_usage(os);
        }

	private:
		void init(const Size &ls);
		bool solve_preconditioned(const Operator<Vector> &A, const Vector &b, Vector &x);
		bool solve_unpreconditioned(const Operator<Vector> &A, const Vector &b, Vector &x);

		Vector r0_;
		Vector r_;
		Vector v_;
		Vector p_;
		Vector h_;
		Vector y_;
		Vector t_;
		Vector s_;
		Vector z_;
		Vector K_inv_t_;
	};
}

#endif //UTOPIA_BICG_STAB_HPP

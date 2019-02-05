#include "utopia_petsc_TaoQPSolver.hpp"
#include "utopia_petsc_TaoSolver.hpp"

#include "utopia_QuadraticFunction.hpp"
#include "utopia_make_unique.hpp"

#include "petsctao.h"

#include <cassert>

namespace utopia {

	template<class Matrix, class Vector>
	class TaoQPSolver<Matrix, Vector>::Impl {
	public:
		typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

		Impl(const std::shared_ptr<LinearSolver> &linear_solver)
		: tao(utopia::make_unique<TaoSolver<Matrix, Vector>>(linear_solver))
		{
			tao->set_type(TAOGPCG);
		}

		std::unique_ptr<TaoSolver<Matrix, Vector>> tao;
	};


	template<class Matrix, class Vector>
	TaoQPSolver<Matrix, Vector>::TaoQPSolver(const std::shared_ptr<LinearSolver> &linear_solver)
	: impl_(utopia::make_unique<Impl>(linear_solver))
	{ }

	template<class Matrix, class Vector>
	TaoQPSolver<Matrix, Vector>::~TaoQPSolver()
	{ }

	template<class Matrix, class Vector>
	TaoQPSolver<Matrix, Vector> * TaoQPSolver<Matrix, Vector>::clone() const {
		m_utopia_warning_once("TaoQPSolver * TaoQPSolver::clone() not implemented properly");
		//FIXME
		//exception-safe cloning
		// auto cloned_tao    = std::unique_ptr<TaoSolver<Matrix, Vector>>(impl_->tao->clone());
		auto cloned        = utopia::make_unique<TaoQPSolver>();
		// cloned->impl_->tao = std::move(cloned_tao);
		return cloned.release();
	}

	template<class Matrix, class Vector>
	bool TaoQPSolver<Matrix, Vector>::apply(const Vector &rhs, Vector &sol)
	{
		assert(this->has_operator());

		impl_->tao->atol(this->atol());
		impl_->tao->stol(this->stol());
		impl_->tao->rtol(this->rtol());
		impl_->tao->max_it(this->max_it());

		impl_->tao->set_box_constraints(this->get_box_constraints());

		QuadraticFunction<Matrix, Vector> fun(
			std::make_shared<Matrix>(*this->get_operator()),
			std::make_shared<Vector>(rhs)
		);

		return impl_->tao->solve(fun, sol);
	}

	template<class Matrix, class Vector>
	void TaoQPSolver<Matrix, Vector>::pc_type(const std::string & pc_type)
	{
	    impl_->tao->set_pc_type(pc_type);
	}

	template<class Matrix, class Vector>
	void TaoQPSolver<Matrix, Vector>::tao_type(const std::string &type)
	{
	    impl_->tao->set_type(type);
	}

	template<class Matrix, class Vector>
	void TaoQPSolver<Matrix, Vector>::read(Input &in)
	{
		impl_->tao->read(in);
	}

	template<class Matrix, class Vector>
	void TaoQPSolver<Matrix, Vector>::set_ksp_types(const std::string &ksp_type, const std::string &pc_type, const std::string &solver_package)
	{
	    impl_->tao->set_ksp_types(ksp_type, pc_type, solver_package);
	}

}

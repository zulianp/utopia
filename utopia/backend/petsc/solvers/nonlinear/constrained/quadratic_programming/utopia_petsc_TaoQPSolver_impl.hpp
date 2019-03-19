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
		auto cloned_tao    = std::unique_ptr<TaoSolver<Matrix, Vector>>(impl_->tao->clone());
		auto cloned        = utopia::make_unique<TaoQPSolver>();
		cloned->impl_->tao = std::move(cloned_tao);
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
	void TaoQPSolver<Matrix, Vector>::set_linear_solver(const std::shared_ptr<LinearSolver> &linear_solver)
    {
        impl_->tao->set_linear_solver(linear_solver);
    }

	template<class Matrix, class Vector>
   	typename TaoQPSolver<Matrix, Vector>::Scalar TaoQPSolver<Matrix, Vector>::atol() const
    {
    	return impl_->tao->atol();
    }

	template<class Matrix, class Vector>
    typename TaoQPSolver<Matrix, Vector>::Scalar TaoQPSolver<Matrix, Vector>::rtol() const
    {
    	return impl_->tao->rtol();
    }

	template<class Matrix, class Vector>
    typename TaoQPSolver<Matrix, Vector>::Scalar TaoQPSolver<Matrix, Vector>::stol() const
    {
    	return impl_->tao->stol();
    }

    template<class Matrix, class Vector>
    typename TaoQPSolver<Matrix, Vector>::SizeType TaoQPSolver<Matrix, Vector>::max_it() const
    {
    	return impl_->tao->max_it();
    }

	template<class Matrix, class Vector>
    bool TaoQPSolver<Matrix, Vector>::verbose() const
    {
    	return impl_->tao->verbose();
    }

	template<class Matrix, class Vector>
    void TaoQPSolver<Matrix, Vector>::atol(const Scalar &atol)
    {
    	impl_->tao->atol(atol);
    }

	template<class Matrix, class Vector>
    void TaoQPSolver<Matrix, Vector>::rtol(const Scalar &rtol)
    {
    	impl_->tao->rtol(rtol);
    }

	template<class Matrix, class Vector>
    void TaoQPSolver<Matrix, Vector>::stol(const Scalar &stol)
    {
    	impl_->tao->stol(stol);
    }

	template<class Matrix, class Vector>
    void TaoQPSolver<Matrix, Vector>::max_it(const SizeType & max_it)
    {
    	impl_->tao->max_it(max_it);
    }

	template<class Matrix, class Vector>
    void TaoQPSolver<Matrix, Vector>::verbose(const bool &verbose)
    {
    	impl_->tao->verbose(verbose);
    }

}

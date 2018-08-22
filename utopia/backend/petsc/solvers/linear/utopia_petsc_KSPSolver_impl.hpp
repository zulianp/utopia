#ifndef UTOPIA_PETSC_KSP_SOLVER_IMPL_HPP
#define UTOPIA_PETSC_KSP_SOLVER_IMPL_HPP

#include "utopia_petsc_KSPSolver.hpp"

#include "utopia_Preconditioner.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_petsc_KSPWrapper.hpp"

#include <algorithm>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsys.h>

namespace utopia {

	class KSPLog {
	public:
		Vec x_k_1;
		Vec x_k_2;

		KSPLog()
		: x_k_1(nullptr), x_k_2(nullptr)
		{}

		inline bool initialized() const
		{
			return x_k_1 != nullptr;
		}

		void init_from(Vec x)
		{
			destroy();

			VecDuplicate(x, &x_k_1);
			VecDuplicate(x, &x_k_2);
		}

		~KSPLog()
		{
			destroy();
		}

		inline void destroy()
		{
			if(x_k_1) {
				VecDestroy(&x_k_1);
				x_k_1 = nullptr;
			}

			if(x_k_2) {
				VecDestroy(&x_k_2);
				x_k_2 = nullptr;
			}
		}
	};

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC>::KSPSolver(const Parameters params)
	: ksp_(std::make_shared<KSPWrapper<Matrix, Vector>>(PETSC_COMM_WORLD, params))
	{
		ksp_type("bcgs");
		pc_type("jacobi");
		ksp_->set_initial_guess_non_zero(true);
		set_parameters(params);
	}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC>::KSPSolver(const std::shared_ptr<KSPWrapper<Matrix, Vector>> &w)
	: ksp_(w)
	{}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::wrap(KSP &ksp)
	{
		ksp_ = std::make_shared< KSPWrapper<Matrix, Vector> >(ksp, false);
	}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC>::~KSPSolver() {}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::set_parameters(const Parameters params)
	{
		PreconditionedSolver::set_parameters(params);
		ksp_->set_parameters(params);
	}

	template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::set_initial_guess_non_zero(const bool val)
	{
		ksp().set_initial_guess_non_zero(val);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::pc_type(const std::string &pc_type)
	{
		ksp_->pc_type(pc_type);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::ksp_type(const std::string & ksp_type)
	{
		ksp_->ksp_type(ksp_type);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::solver_package(const std::string &package)
	{
		ksp_->solver_package(package);
	}

    template<typename Matrix, typename Vector>
	std::string KSPSolver<Matrix, Vector, PETSC>::pc_type() const { return this->ksp_->pc_type();}

    template<typename Matrix, typename Vector>
	std::string KSPSolver<Matrix, Vector, PETSC>::ksp_type() const { return this->ksp_->ksp_type();}

    template<typename Matrix, typename Vector>
	std::string KSPSolver<Matrix, Vector, PETSC>::solver_package() const { return this->ksp_->solver_package();}

    template<typename Matrix, typename Vector>
	bool KSPSolver<Matrix, Vector, PETSC>::apply(const Vector &b, Vector &x)
	{
		ksp_->set_tolerances(this->rtol(), this->atol(), PETSC_DEFAULT, this->max_it());

        // is this proper place to do so???
        // this->set_ksp_options(ksp_->implementation());
		return ksp_->apply(b, x);
	}

    template<typename Matrix, typename Vector>
	bool KSPSolver<Matrix, Vector, PETSC>::smooth(const Vector &rhs, Vector &x)
	{
		return ksp_->smooth(this->sweeps(), rhs, x);
	}

    template<typename Matrix, typename Vector>
	bool KSPSolver<Matrix, Vector, PETSC>::must_compute_cond_number() const
	{
		static const std::vector<std::string> types = {"bcgs", "cg", "gmres"};
		return (this->verbose() && in_array(ksp_->ksp_type(), types));
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::set_ksp_options(KSP &ksp)
	{
		ksp_->copy_settings_to(ksp);
		set_monitor_options(ksp);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::attach_preconditioner(KSP &ksp) const {
		KSPWrapper<Matrix, Vector> w(ksp, false);
		auto delegate_ptr = std::dynamic_pointer_cast<DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());

		if(delegate_ptr) {
			if(ksp_->has_shell_pc()) {
				m_utopia_warning_once("set_preconditioner sets jacobi if a delegate precond has been set and type is matshell");
			}
		} else if(this->get_preconditioner()) {
			auto shell_ptr = this->get_preconditioner().get();
			w.attach_shell_preconditioner(UtopiaPCApplyShell,
				shell_ptr,
				nullptr,
				nullptr
				);
		}
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
	{
		PreconditionedSolver::set_preconditioner(precond);

		auto delegate_ptr = std::dynamic_pointer_cast<DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());

		if(delegate_ptr) {
			if(ksp_->has_shell_pc()) {
				m_utopia_warning_once("set_preconditioner sets jacobi if a delegate precond has been set and type is matshell");
				ksp_->pc_type("jacobi");
			}

		} else if(this->get_preconditioner()) {
			auto shell_ptr = this->get_preconditioner().get();
			ksp_->attach_shell_preconditioner(UtopiaPCApplyShell,
				shell_ptr,
				nullptr,
				nullptr
				);
		}
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::set_monitor_options(KSP &ksp) const
	{
		PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
		if(this->verbose()) {
			ierr = KSPMonitorSet(ksp,
				MyKSPMonitor,
				new KSPLog(),
				[](void **arg) -> PetscErrorCode {
					auto ksp_log = (KSPLog **)(arg);
					delete *ksp_log;
					*ksp_log = nullptr;
					return 0;
				});  assert(ierr == 0);
		}

		if(must_compute_cond_number()) {
			ierr = KSPSetComputeSingularValues(ksp, PETSC_TRUE); assert(ierr == 0);
		}
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::handle_reset(const Matrix &op)
	{
		bool must_reset = true;
        // bool must_reset = false;

        // if(this->has_operator()) {
        //     auto s = size(*this->get_operator());
        //     auto other_s = size(op);
        //     if(s != other_s) {
        //     must_reset = true;
        //     }
        // }

        // if(op.implementation().communicator() != ksp_->communicator())
        // {
        //     must_reset = true;
        // }

		if(must_reset) {
			auto temp_ksp = std::make_shared<KSPWrapper<Matrix, Vector>>(op.implementation().communicator());
			temp_ksp->copy_settings_from(*ksp_);
			ksp_ = temp_ksp;

			if(this->get_preconditioner()) {
				set_preconditioner(this->get_preconditioner());
			}
		}
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec)
	{
		handle_reset(*op);
		set_monitor_options(ksp_->implementation());


		PreconditionedSolver::update(op, prec);
		ksp_->update(*op, *prec);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::update(const std::shared_ptr<const Matrix> &op)
	{
		handle_reset(*op);

		PreconditionedSolver::update(op);
		set_monitor_options(ksp_->implementation());

		bool skip_set_operators = false;
		if(this->get_preconditioner()) {
			auto mat_prec = std::dynamic_pointer_cast< DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());
			if(mat_prec) {
				ksp_->update(*op, *mat_prec->get_matrix());
				skip_set_operators = true;
			}
		}

		if(!skip_set_operators) {
			ksp_->update(*op);
		}
	}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC> & KSPSolver<Matrix, Vector, PETSC>::operator=(const KSPSolver<Matrix, Vector, PETSC> &other)
	{
		if(this == &other) return *this;
		PreconditionedSolver::operator=(other);
		Smoother::operator=(other);

		ksp_ = std::make_shared<KSPWrapper<Matrix, Vector>>(other.ksp_->communicator());
		ksp_->copy_settings_from(*other.ksp_);
		return *this;
	}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC> & KSPSolver<Matrix, Vector, PETSC>::operator=(KSPSolver<Matrix, Vector, PETSC> &&other)
	{
		if(this == &other) return *this;
		PreconditionedSolver::operator=(std::move(other));
		Smoother::operator=(std::move(other));
		ksp_ = std::move(other.ksp_);
		return *this;
	}


    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC>::KSPSolver(const KSPSolver<Matrix, Vector, PETSC> &other):
	PreconditionedSolver(other),
	Smoother(other),
	ksp_(std::make_shared<KSPWrapper<Matrix, Vector>>(other.ksp_->communicator()))
	{
		ksp_->copy_settings_from(*other.ksp_);
	}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC>::KSPSolver(KSPSolver<Matrix, Vector, PETSC> &&other)
	: PreconditionedSolver(std::move(other)),
	Smoother(std::move(other)),
	ksp_(std::move(other.ksp_))
	{}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC> * KSPSolver<Matrix, Vector, PETSC>::clone() const
	{
		return new KSPSolver(*this);
	}

    template<typename Matrix, typename Vector>
	KSPWrapper<Matrix, Vector> &KSPSolver<Matrix, Vector, PETSC>::ksp()
	{
		return *ksp_;
	}

    template<typename Matrix, typename Vector>
	const KSPWrapper<Matrix, Vector> &KSPSolver<Matrix, Vector, PETSC>::ksp() const
	{
		return *ksp_;
	}

	template<typename Matrix, typename Vector>
	KSP &KSPSolver<Matrix, Vector, PETSC>::implementation()
	{
		return ksp().implementation();
	}
}


#endif //UTOPIA_PETSC_KSP_SOLVER_IMPL_HPP

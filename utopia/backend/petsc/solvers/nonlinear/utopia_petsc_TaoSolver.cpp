#include "utopia_petsc_TaoSolver.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Types.hpp"
#include "utopia_Describable.hpp"

#include "petsctao.h"
#include <mpi.h>

#define U_CHECKERR(ierr) { if(ierr != 0) return false; }

namespace utopia {

    class TaoTypes : public Describable {
    public:

        inline static bool is_valid(const std::string &type, const bool verbose = true)
        {
            const auto &i = instance();
            bool valid = i.types_.find(type) != i.types_.end();

            if(!valid && verbose) {
                std::cerr << "Invalid tao type " << type << ". Valid types are: " << std::endl;
                i.describe(std::cerr);
            }

            return valid;
        }

        void describe(std::ostream &os) const override
        {
            for(const auto &t : types_) {
                os << t << " ";
            }

            os << std::endl;
        }


    private:
        std::set<std::string> types_;

        static inline const TaoTypes &instance()
        {
            static TaoTypes instance_;
            return instance_;
        }

        TaoTypes()
        {
            types_.insert(TAOLMVM);
            types_.insert(TAONLS);
            types_.insert(TAONTR);
            types_.insert(TAONTL);
            types_.insert(TAOCG);
            types_.insert(TAOTRON);
            types_.insert(TAOOWLQN);
            types_.insert(TAOBMRM);
            types_.insert(TAOBLMVM);

            types_.insert(TAOBQPIP);
            types_.insert(TAOGPCG);
            types_.insert(TAONM);
            types_.insert(TAOPOUNDERS);
            types_.insert(TAOLCL);
            types_.insert(TAOSSILS);
            types_.insert(TAOSSFLS);
            types_.insert(TAOASILS);
            types_.insert(TAOASFLS);
            types_.insert(TAOIPM);

            //TODO check if this is the right version
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 10, 3)
          
            types_.insert(TAOBQNLS);
            types_.insert(TAOBNCG);

            types_.insert(TAOBNLS);
            types_.insert(TAOBNTR);
            types_.insert(TAOBNTL);

            types_.insert(TAOBQNKLS);
            types_.insert(TAOBQNKTR);
            types_.insert(TAOBQNKTL);
#endif
        }
    };

    
    template<class Matrix, class Vector>
    static PetscErrorCode UtopiaTaoEvaluateObjective(Tao tao, Vec x, PetscReal *ret, void *ctx)
    {
    	UTOPIA_UNUSED(tao);
        using Scalar = typename Function<Matrix, Vector>::Scalar;
        
        Function<Matrix, Vector> * fun = static_cast<Function<Matrix, Vector> *>(ctx);
        assert(fun);
        
        Vector utopia_x;
        convert(x, utopia_x);
        
        Scalar utopia_value = 0.;
        if(!fun->value(utopia_x, utopia_value)) {
            return 1;
        }
        
        *ret = utopia_value;
        return 0;
    }
    
    template<class Matrix, class Vector>
    static PetscErrorCode UtopiaTaoEvaluateGradient(Tao tao, Vec x, Vec g, void *ctx)
    {
        UTOPIA_UNUSED(tao);
        Function<Matrix, Vector> * fun = static_cast<Function<Matrix, Vector> *>(ctx);
        assert(fun);
        
        Vector utopia_x;
        convert(x, utopia_x);
        
        Vector utopia_g = local_zeros(local_size(utopia_x));
        if(!fun->gradient(utopia_x, utopia_g)) {
            return 1;
        }
        
        convert(utopia_g, g);
        return 0;
    }
    
    template<class Matrix, class Vector>
    static PetscErrorCode UtopiaTaoFormHessian(Tao tao, Vec x, Mat H, Mat Hpre, void *ctx)
    {
        UTOPIA_UNUSED(tao);
        // PetscErrorCode ierr = 0;
        Function<Matrix, Vector> * fun = static_cast<Function<Matrix, Vector> *>(ctx);
        assert(fun);
        
        Vector utopia_x;
        
        Matrix utopia_H;
        Matrix utopia_Hpre;
        
        convert(x, utopia_x);
        utopia_H.implementation().wrap(H);
        utopia_Hpre.implementation().wrap(Hpre);
        
        if(!fun->hessian(utopia_x, utopia_H, utopia_Hpre)) {
            if(!fun->hessian(utopia_x, utopia_H)) {
                utopia_error("[Error] Failed to assemble Hessian.");
                assert(false);
                return 1;
            }
            
        } else {
            if(raw_type(utopia_Hpre) != Hpre) {
                //FIXME maybe add optimization options
                MatCopy(raw_type(utopia_Hpre), Hpre, DIFFERENT_NONZERO_PATTERN);
            }
        }
        
        if(raw_type(utopia_H) != H) {
            //FIXME maybe add optimization options
            MatCopy(raw_type(utopia_H), H, DIFFERENT_NONZERO_PATTERN);
        }
        
        return 0;
    }
    
    template<class Matrix, class Vector>
    bool UtopiaTaoSetUp(Tao tao, Function<Matrix, Vector> &fun)
    {
        fun.data()->init();
        if(!fun.initialize_hessian(*fun.data()->H, *fun.data()->H_pre)) {
            utopia_error("TaoSolver requires Function::initialize_hessian to be implemented.");
            assert(false);
            return false;
        }
        
        PetscErrorCode ierr = 0;
        if(fun.has_preconditioner())
        {
            ierr = TaoSetHessianRoutine(tao,
                                        raw_type(*fun.data()->H),
                                        raw_type(*fun.data()->H_pre),
                                        UtopiaTaoFormHessian<Matrix, Vector>,
                                        static_cast<void *>(&fun)); U_CHECKERR(ierr);
            
        } else {
            ierr = TaoSetHessianRoutine(tao,
                                        raw_type(*fun.data()->H),
                                        raw_type(*fun.data()->H),
                                        UtopiaTaoFormHessian<Matrix, Vector>,
                                        static_cast<void *>(&fun)); U_CHECKERR(ierr);
        }
        
        ierr = TaoSetObjectiveRoutine(tao,
                                      UtopiaTaoEvaluateObjective<Matrix, Vector>,
                                      static_cast<void *>(&fun)); U_CHECKERR(ierr);
        
        ierr = TaoSetGradientRoutine(
                                     tao,
                                     UtopiaTaoEvaluateGradient<Matrix, Vector>,
                                     static_cast<void *>(&fun)); U_CHECKERR(ierr);
        
        return false;
    }
    
    template<class Matrix, class Vector>
    class TaoSolver<Matrix, Vector>::Impl : public Configurable {
    public:	

    	using Scalar   = UTOPIA_SCALAR(Vector);
    	using SizeType = UTOPIA_SIZE_TYPE(Vector);

    	Impl(MPI_Comm comm)
    	: tao(nullptr), type_(TAOTRON)
    	{
            assert(TaoTypes::is_valid(type_));
    		init(comm);
    	}

    	Impl()
    	: tao(nullptr), type_(TAOTRON)
    	{
            assert(TaoTypes::is_valid(type_));
        }

    	~Impl()
    	{}

    	void init(MPI_Comm comm)
    	{
    		destroy();

            assert(TaoTypes::is_valid(type_));

    		auto ierr = TaoCreate(comm, &tao);  assert(ierr == 0);
    		ierr      = TaoSetType(tao, type_.c_str());	assert(ierr == 0);
    	}

    	void set_from_options()
    	{
    		assert(initialized());    
    		auto ierr = TaoSetFromOptions(tao);     assert(ierr == 0);
    	}

    	void destroy()
    	{
    		if(tao) {
    			TaoDestroy(&tao);
    			tao = nullptr;
    		}
    	}

    	inline std::string get_type() const
    	{
            if(initialized()) {
                //TODO check if this is the right version
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 10, 2)
                TaoType type;
#else
                const TaoType type;
#endif
                TaoGetType(tao, &type);

                assert(TaoTypes::is_valid(type));
                return type;
            }

            assert(TaoTypes::is_valid(type_));

    		return type_;
    	}

    	void set_type(const std::string &type)
    	{
            assert(TaoTypes::is_valid(type));

    		type_ = type;
    		if(tao) {
    			TaoSetType(tao, type_.c_str());
    		}
    	}

    	void set_bounds(const Vector &lb, const Vector &ub)
    	{    	
    		assert(initialized());    

    		if(!tao) {
    			std::cerr << "[Error] attempt to set bounds to uninitialized tao" << std::endl;
    			return;
    		}

    	    auto ierr = TaoSetVariableBounds(tao, raw_type(lb), raw_type(ub)); assert(ierr == 0);
    	}

    	bool get_ksp(KSP *ksp)
	    {
	        assert(initialized());

	        if(!tao) {
	        	std::cerr << "[Error] attempt to set ksp to uninitialized tao" << std::endl;
	        	return false;
	        }
	        
	        PetscErrorCode ierr = TaoGetKSP(tao, ksp); assert(ierr == 0);
	        return ierr == 0;
	    }

    	void set_tol(const Scalar gatol,
                     const Scalar grtol,
					 const Scalar gttol,
                     const SizeType maxits)
    	{
    		assert(initialized());

    		if(!tao) {
    			std::cerr << "[Error] attempt to set tol to uninitialized tao" << std::endl;
    			return;
    		}

    		auto ierr = TaoSetTolerances(tao, gatol, grtol, gttol); assert(ierr == 0);
    		ierr = TaoSetMaximumIterations(tao, maxits); 			assert(ierr == 0);
    	}

    	void set_monitor()
	    {
	        const char  monfilename[7] ="stdout";
	        PetscViewer    monviewer;
	        PetscViewerASCIIOpen(communicator(), monfilename, &monviewer);
	        TaoSetMonitor(tao, TaoDefaultSMonitor, monviewer, (PetscErrorCode (*)(void**))PetscViewerDestroy);
	    }

	    void read(Input &in) override
	    {
	    	std::string type;
	    	in.get("type", type);
	    	set_type(type.c_str());
	    }

	    inline bool initialized() const
	    {
	    	return tao != nullptr;
	    }

	    void set_linear_solver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &solver)
	    {
	    	assert(tao);

	    	KSP ksp;
	    	get_ksp(&ksp);
	    	build_ksp(solver, ksp);
	    }

	   void set_function(Function<Matrix, Vector> &fun)
        {
            UtopiaTaoSetUp(tao, fun);
        }


	    inline bool solve(Vector &x)
	    {
    	    PetscErrorCode ierr = 0;
    	    TaoSetInitialVector(tao, raw_type(x));
    	    ierr = TaoSolve(tao); U_CHECKERR(ierr);
    	    
    	    PetscInt iterate;
    	    PetscReal f;
    	    PetscReal gnorm;
    	    PetscReal cnorm;
    	    PetscReal xdiff;
    	    TaoConvergedReason reason;
    	    TaoGetSolutionStatus(tao, &iterate, &f, &gnorm, &cnorm, &xdiff, &reason);
    	    
    	    // if(this->verbose()) {
    	    // std::cout << "iterate: " << iterate << std::endl;
    	    // std::cout << "f: " << f << std::endl;
    	    // std::cout << "gnorm: " << gnorm << std::endl;
    	    // std::cout << "cnorm: " << cnorm << std::endl;
    	    // std::cout << "xdiff: " << xdiff << std::endl;
    	    // std::cout << "reason: " << reason << std::endl;
    	    // }
    	    
    	    if(reason < 0) {
    	        utopia_warning("> Failed to converge");
    	    }
    	    
    	    return reason >= 0;
	    }

	    inline bool smooth(Vector &x)
	    {
	    	PetscErrorCode ierr = 0;
	    	TaoSetInitialVector(tao, raw_type(x));
	    	ierr = TaoSolve(tao); U_CHECKERR(ierr);
	    	return true;
	    }
    

	    inline MPI_Comm communicator() const {
	    	assert(initialized());
	    	MPI_Comm comm = PetscObjectComm((PetscObject) tao);
	    	assert(comm != MPI_COMM_NULL);
	    	return comm;
	    }

    private:
    	Tao tao;
    	std::string type_;

    };
     
    template<class Matrix, class Vector>
    TaoSolver<Matrix, Vector>::TaoSolver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &linear_solver)
    : NewtonBase<Matrix, Vector>(linear_solver)
    {
    	impl_ = utopia::make_unique<Impl>();
    }
    
    template<class Matrix, class Vector>
    TaoSolver<Matrix, Vector>::TaoSolver()
    : NewtonBase<Matrix, Vector>(nullptr)
    {
    	impl_ = utopia::make_unique<Impl>();
    }
    
    template<class Matrix, class Vector>
    TaoSolver<Matrix, Vector>::~TaoSolver() {}
    
    template<class Matrix, class Vector>
    void TaoSolver<Matrix, Vector>::set_type(const std::string &type)
    {
       impl_->set_type(type.c_str());
    }
    
    template<class Matrix, class Vector>
    void TaoSolver<Matrix, Vector>::read(Input &in)
    {
        NewtonBase<Matrix, Vector>::read(in);
        // VariableBoundSolverInterface<Vector>::read(in);
     	impl_->read(in);
    }
    
    
    template<class Matrix, class Vector>
    void TaoSolver<Matrix, Vector>::print_usage(std::ostream &os) const
    {
        NewtonBase<Matrix, Vector>::print_usage(os);
        this->print_param_usage(os, "type", "string", "Type of tao solver.", "-");
    }
    
    template<class Matrix, class Vector>
    bool TaoSolver<Matrix, Vector>::solve(Function<Matrix, Vector> &fun, Vector &x)
    {
        init(fun, x);
        return impl_->solve(x);
    }
    
    template<class Matrix, class Vector>
    bool TaoSolver<Matrix, Vector>::smooth(Function<Matrix, Vector> &fun, Vector &x)
    {
        init(fun, x);
        return impl_->smooth(x);
    }
    
    template<class Matrix, class Vector>
    void TaoSolver<Matrix, Vector>::init(Function<Matrix, Vector> &fun, Vector & x)
    {
    	if(!impl_->initialized()) {
    		impl_->init(x.implementation().communicator());

    		 if(this->linear_solver()) {
        		impl_->set_linear_solver(this->linear_solver());
        	}
    	}
        
        impl_->set_tol(this->atol(), this->rtol(), this->stol(), this->max_it());

        if(this->has_bound()) {
            this->fill_empty_bounds();
            impl_->set_bounds(this->get_lower_bound(), this->get_upper_bound());
        }

        impl_->set_function(fun);

        if(this->verbose()) {
            impl_->set_monitor();
        }
        
        impl_->set_from_options();
    }
    
    template<class Matrix, class Vector>
    bool TaoSolver<Matrix, Vector>::get_ksp(KSP *ksp)
    {
        return impl_->get_ksp(ksp);
    }

	template<class Matrix, class Vector>
    TaoSolver<Matrix, Vector> *TaoSolver<Matrix, Vector>::clone() const
    {
    	//FIXME make proper deep clone
    	auto cloned = utopia::make_unique<TaoSolver>();

    	cloned->set_type(cloned->impl_->get_type());
    	
    	if(this->impl_->initialized()) {
    		cloned->impl_->init(impl_->communicator());
    	}

    	return cloned.release();
    }
    
    template class TaoSolver<DSMatrixd, DVectord>;
    template class TaoSolver<DMatrixd, DVectord>;
}

#undef U_CHECKERR

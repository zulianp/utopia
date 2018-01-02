#ifndef UTOPIA_PETSC_EVAL_FACTORY_HPP
#define UTOPIA_PETSC_EVAL_FACTORY_HPP 

namespace utopia {
	template<class Index, class Traits>
	class Eval< Ghosts<Index>, Traits, PETSC> {
	public:
		typedef typename TypeAndFill<Traits, Ghosts<Index> >::Type Return;

	    inline static Return apply(const Ghosts<Index> &expr) {
	        Return ret;

	        UTOPIA_LOG_BEGIN(expr);

	        UTOPIA_BACKEND(Traits).build_ghosts(expr.local_size(), expr.global_size(), expr.index(), ret);

	        UTOPIA_LOG_END(expr);
	        return ret;
	    }
	};
}

#endif //UTOPIA_PETSC_EVAL_FACTORY_HPP

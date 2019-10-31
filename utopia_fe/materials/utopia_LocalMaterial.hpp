#ifndef UTOPIA_LOCAL_MATERIAL_HPP
#define UTOPIA_LOCAL_MATERIAL_HPP

#include "utopia_fe_base.hpp"

namespace utopia {
	
	template<class Fun>
	inline void loop(const SizeType n, Fun f)
	{
	    for(SizeType i = 0; i < n; ++i) {
	        f(i);
	    }
	}

	template<class FE>
	class LocalMaterial {
	public:
		virtual ~LocalMaterial() {}

	    virtual void init(FE &element) = 0;
	   
	    virtual void assemble(FE &element, USerialMatrix &mat) = 0;
	    virtual void assemble(FE &element, USerialVector &vec) = 0;

	    virtual void assemble(FE &element, USerialMatrix &mat, USerialVector &vec)
	    {
	        assemble(element, mat);
	        assemble(element, vec);
	    }
	};
}

#endif

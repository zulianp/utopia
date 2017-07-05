#ifndef UTOPIA_S_TREE_HPP
#define UTOPIA_S_TREE_HPP 

#include "utopia_SElementAdapter.hpp"

#include "SparseExpression.hpp"
#include "express_Sparse.hpp"

#include "cutlibpp.hpp"
#include "cutlibpp_Base.hpp"
#include "cutlibpp_Tree.hpp"
#include "cutlibpp_NTreeMutatorFactory.hpp"
#include "cutlibpp_NTreeWithSpanMutatorFactory.hpp"
#include "cutlibpp_NTreeWithTagsMutatorFactory.hpp"
#include "cutlibpp_API.hpp"

namespace utopia {
	template<int _Dimension>
	class STreeTraits {
	public:
	    enum {
	        Dimension = _Dimension
	    };
	    
	    typedef utopia::BoxBoxAdapter<Dimension> Bound;
	    typedef utopia::SElementAdapter<Dimension> DataType;  
	};
	    
	template<int Dimension>
	class STree : public cutlibpp::Tree< STreeTraits<Dimension> > {
	public:
	    typedef STreeTraits<Dimension> Traits;
	    
	    STree() {};
	    
	    static cutk::shared_ptr<STree> New(const cutk::shared_ptr<cutlibpp::Predicate> &predicate,
	                                                    const int maxElementsXNode = cutlibpp::DEFAULT_REFINE_MAX_ELEMENTS,
	                                                    const int maxDepth = cutlibpp::DEFAULT_REFINE_DEPTH
	                                                    ) {
	        using namespace cutlibpp;
	        cutk::shared_ptr<STree> tree = cutk::make_shared<STree>();
	        cutk::shared_ptr<NTreeWithTagsMutatorFactory < STree> > factory =
	        cutk::make_shared<NTreeWithTagsMutatorFactory < STree> >(predicate);
	        factory->setRefineParams(maxElementsXNode, maxDepth);
	        tree->setMutatorFactory(factory);
	        return tree;
	    } 
	};
}


#endif //UTOPIA_S_TREE_HPP


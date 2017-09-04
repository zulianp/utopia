#ifndef UTOPIA_S_TREE_HPP
#define UTOPIA_S_TREE_HPP 

#include "utopia_SElementAdapter.hpp"

#include "moonolith_tree.hpp"
#include "moonolith_n_tree_mutator_factory.hpp"
#include "moonolith_n_tree_with_span_mutator_factory.hpp"
#include "moonolith_n_tree_with_tags_mutator_factory.hpp"
#include "par_moonolith.hpp"

#include <memory>

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
	class STree : public moonolith::Tree< STreeTraits<Dimension> > {
	public:
	    typedef STreeTraits<Dimension> Traits;
	    
	    STree() {};
	    
	    static std::shared_ptr<STree> New(const std::shared_ptr<moonolith::Predicate> &predicate,
	                                                    const int maxElementsXNode = moonolith::DEFAULT_REFINE_MAX_ELEMENTS,
	                                                    const int maxDepth = moonolith::DEFAULT_REFINE_DEPTH
	                                                    ) {
	        using namespace moonolith;
	        std::shared_ptr<STree> tree = std::make_shared<STree>();
	        std::shared_ptr<NTreeWithTagsMutatorFactory<STree>> factory =
	        std::make_shared<NTreeWithTagsMutatorFactory<STree>>(predicate);
	        factory->set_refine_params(maxElementsXNode, maxDepth);
	        tree->set_mutator_factory(factory);
	        return tree;
	    } 
	};
}


#endif //UTOPIA_S_TREE_HPP


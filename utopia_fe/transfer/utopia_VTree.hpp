
#ifndef UTOPIA_V_TREE_HPP
#define UTOPIA_V_TREE_HPP 

#include "utopia_SElementAdapter.hpp"

#include "moonolith_tree.hpp"
#include "moonolith_n_tree_mutator_factory.hpp"
#include "moonolith_n_tree_with_span_mutator_factory.hpp"
#include "moonolith_n_tree_with_tags_mutator_factory.hpp"
#include "moonolith_api.hpp"

namespace utopia {
 template<int _Dimension>
    class VTreeTraits {
    public:
        enum {
            Dimension = _Dimension
        };
        
        typedef utopia::BoxBoxAdapter<Dimension> Bound;
        typedef utopia::VElementAdapter<Dimension> DataType;
        
    };
    
    template<int Dimension>
    class VTree : public moonolith::Tree< VTreeTraits<Dimension> > {
    public:
        typedef VTreeTraits<Dimension> Traits;
        
        VTree() {};
        
        
        static std::shared_ptr<VTree> New(const std::shared_ptr<moonolith::Predicate> &predicate,
                                                 const int maxElementsXNode=moonolith::DEFAULT_REFINE_MAX_ELEMENTS,
                                                 const int maxDepth=moonolith::DEFAULT_REFINE_DEPTH
                                                 ) {
            using namespace cutlibpp;
            std::shared_ptr<VTree> tree = moonolith::make_shared<VTree>();
            std::shared_ptr<NTreeWithTagsMutatorFactory < VTree> > factory =
            moonolith::make_shared<NTreeWithTagsMutatorFactory < VTree> >(predicate);
            factory->setRefineParams(maxElementsXNode, maxDepth);
            std::cout<<"LibMeshTree::maxElements = "<<maxElementsXNode<<std::endl;
            std::cout<<"LibMeshTree::maxDepth = "<<maxDepth<<std::endl;
            tree->setMutatorFactory(factory);
            return tree;
        }
        
        
    };
}


#endif //UTOPIA_V_TREE_HPP


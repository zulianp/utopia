
#ifndef UTOPIA_V_TREE_HPP
#define UTOPIA_V_TREE_HPP 

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
    class VTreeTraits {
    public:
        enum {
            Dimension = _Dimension
        };
        
        typedef utopia::BoxBoxAdapter<Dimension> Bound;
        typedef utopia::VElementAdapter<Dimension> DataType;
        
    };
    
    template<int Dimension>
    class VTree : public cutlibpp::Tree< VTreeTraits<Dimension> > {
    public:
        typedef VTreeTraits<Dimension> Traits;
        
        VTree() {};
        
        
        static cutk::shared_ptr<VTree> New(const cutk::shared_ptr<cutlibpp::Predicate> &predicate,
                                                 const int maxElementsXNode=cutlibpp::DEFAULT_REFINE_MAX_ELEMENTS,
                                                 const int maxDepth=cutlibpp::DEFAULT_REFINE_DEPTH
                                                 ) {
            using namespace cutlibpp;
            cutk::shared_ptr<VTree> tree = cutk::make_shared<VTree>();
            cutk::shared_ptr<NTreeWithTagsMutatorFactory < VTree> > factory =
            cutk::make_shared<NTreeWithTagsMutatorFactory < VTree> >(predicate);
            factory->setRefineParams(maxElementsXNode, maxDepth);
            std::cout<<"LibMeshTree::maxElements = "<<maxElementsXNode<<std::endl;
            std::cout<<"LibMeshTree::maxDepth = "<<maxDepth<<std::endl;
            tree->setMutatorFactory(factory);
            return tree;
        }
        
        
    };
}


#endif //UTOPIA_V_TREE_HPP


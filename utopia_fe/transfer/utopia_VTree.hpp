
#ifndef UTOPIA_V_TREE_HPP
#define UTOPIA_V_TREE_HPP

#include "utopia_VElementAdapter.hpp"

#include "moonolith_tree.hpp"
#include "moonolith_n_tree_mutator_factory.hpp"
#include "moonolith_n_tree_with_span_mutator_factory.hpp"
#include "moonolith_n_tree_with_tags_mutator_factory.hpp"
#include "par_moonolith.hpp"

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
                                          const int maxElementsXNode = moonolith::DEFAULT_REFINE_MAX_ELEMENTS,
                                          const int maxDepth = moonolith::DEFAULT_REFINE_DEPTH
                                                 ) {
            using namespace moonolith;
            auto tree = std::make_shared<VTree>();
            auto factory = std::make_shared<NTreeWithTagsMutatorFactory < VTree> >(predicate);
            factory->set_refine_params(maxElementsXNode, maxDepth);
            tree->set_mutator_factory(factory);
            return tree;
        }


    };
}


#endif //UTOPIA_V_TREE_HPP


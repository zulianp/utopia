#ifndef UTOPIA_TREE_PROPERTIES_HPP
#define UTOPIA_TREE_PROPERTIES_HPP

#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

    constexpr int static_max(const int left, const int right)
    {
        return (left > right)? left : right;
    }

    constexpr int static_min(const int left, const int right)
    {
        return (left < right)? left : right;
    }


    template<class Expr>
    class TreeProperties {
    public:
        enum { greatest_tensor_order = -1 	   };
        enum { smallest_tensor_order = 1000    };
        enum { has_mat_mat_mul 		 = 0 	   };
        enum { has_differentiable_sub_tree = 0 };
    };

    template<class InnerExpr, class Operation>
    class TreeProperties< Unary<InnerExpr, Operation> > {
    public:
        enum { greatest_tensor_order = TreeProperties<InnerExpr>::greatest_tensor_order };
        enum { smallest_tensor_order = TreeProperties<InnerExpr>::smallest_tensor_order };
        enum { has_mat_mat_mul 		 = TreeProperties<InnerExpr>::has_mat_mat_mul 		};
        enum { has_differentiable_sub_tree = TreeProperties<InnerExpr>::has_differentiable_sub_tree };

    };

    // template<class InnerExpr>
    // class TreeProperties< Diag<InnerExpr> > {
    // public:
    // 	enum { greatest_tensor_order = InnerExpr::Order == 2?
    // 									static_max(InnerExpr::Order - 1, TreeProperties<InnerExpr>::greatest_tensor_order) :
    // 									static_max(InnerExpr::Order + 1, TreeProperties<InnerExpr>::greatest_tensor_order) };

    // 	enum { smallest_tensor_order = InnerExpr::Order == 2?
    // 									static_min(InnerExpr::Order - 1, TreeProperties<InnerExpr>::smallest_tensor_order) :
    // 									static_min(InnerExpr::Order + 1, TreeProperties<InnerExpr>::smallest_tensor_order) };

    // 	enum { has_mat_mat_mul 		 = TreeProperties<InnerExpr>::has_mat_mat_mul 		};
    // 	enum { has_differentiable_sub_tree = TreeProperties<InnerExpr>::has_differentiable_sub_tree };

    // };

    template<class InnerExpr, class Operation>
    class TreeProperties< Reduce<InnerExpr, Operation> > {
    public:
        enum { greatest_tensor_order = TreeProperties<InnerExpr>::greatest_tensor_order };
        enum { smallest_tensor_order = 0 												};
        enum { has_mat_mat_mul 		 = TreeProperties<InnerExpr>::has_mat_mat_mul 		};
        enum { has_differentiable_sub_tree = TreeProperties<InnerExpr>::has_differentiable_sub_tree };
    };

    template<class InnerExpr>
    class TreeProperties< Trace<InnerExpr> > {
    public:
        enum { greatest_tensor_order = TreeProperties<InnerExpr>::greatest_tensor_order };
        enum { smallest_tensor_order = 0 												};
        enum { has_mat_mat_mul 		 = TreeProperties<InnerExpr>::has_mat_mat_mul 		};
        enum { has_differentiable_sub_tree = TreeProperties<InnerExpr>::has_differentiable_sub_tree };
    };

    template<class Left, class Right, class Operation>
    class TreeProperties< Binary<Left, Right, Operation> > {
    public:
        enum { greatest_tensor_order = static_max(TreeProperties<Left>::greatest_tensor_order, TreeProperties<Right>::greatest_tensor_order)  };
        enum { smallest_tensor_order = static_min(TreeProperties<Left>::smallest_tensor_order, TreeProperties<Right>::smallest_tensor_order) };
        enum { has_mat_mat_mul = TreeProperties<Left>::has_mat_mat_mul || TreeProperties<Right>::has_mat_mat_mul };
        enum { has_differentiable_sub_tree = TreeProperties<Left>::has_differentiable_sub_tree || TreeProperties<Right>::has_differentiable_sub_tree  };
    };

    template<class Left, class Right>
    class TreeProperties< Multiply<Left, Right> > {
    public:
        enum { greatest_tensor_order = static_max(TreeProperties<Left>::greatest_tensor_order, TreeProperties<Right>::greatest_tensor_order)  };
        enum { smallest_tensor_order = static_min(TreeProperties<Left>::smallest_tensor_order, TreeProperties<Right>::smallest_tensor_order) };
        enum { has_mat_mat_mul = TreeProperties<Left>::has_mat_mat_mul || TreeProperties<Right>::has_mat_mat_mul || ( Left::Order == 2 && Right::Order == 2 ) };
        enum { has_differentiable_sub_tree = TreeProperties<Left>::has_differentiable_sub_tree || TreeProperties<Right>::has_differentiable_sub_tree  };
    };

    template<class Left, class Right>
    class TreeProperties< Construct<Left, Right> > {
    public:
        enum { greatest_tensor_order = static_max(TreeProperties<Left>::greatest_tensor_order, TreeProperties<Right>::greatest_tensor_order)  };
        enum { smallest_tensor_order = static_min(TreeProperties<Left>::smallest_tensor_order, TreeProperties<Right>::smallest_tensor_order) };
        enum { has_mat_mat_mul = TreeProperties<Right>::has_mat_mat_mul };
        enum { has_differentiable_sub_tree = TreeProperties<Right>::has_differentiable_sub_tree  };
    };

    template<class Left, class Right>
    class TreeProperties< Assign<Left, Right> > {
    public:
        enum { greatest_tensor_order = static_max(TreeProperties<Left>::greatest_tensor_order, TreeProperties<Right>::greatest_tensor_order)  };
        enum { smallest_tensor_order = static_min(TreeProperties<Left>::smallest_tensor_order, TreeProperties<Right>::smallest_tensor_order) };
        enum { has_mat_mat_mul = TreeProperties<Left>::has_mat_mat_mul || TreeProperties<Right>::has_mat_mat_mul || ( Left::Order == 2 && Right::Order == 2 ) };
        enum { has_differentiable_sub_tree = TreeProperties<Right>::has_differentiable_sub_tree  };
    };

    template<class InnerExpr>
    class TreeProperties< Transposed<InnerExpr> > {
    public:
        enum { greatest_tensor_order = TreeProperties<InnerExpr>::greatest_tensor_order };
        enum { smallest_tensor_order = TreeProperties<InnerExpr>::smallest_tensor_order };
        enum { has_mat_mat_mul = TreeProperties<InnerExpr>::has_mat_mat_mul };
        enum { has_differentiable_sub_tree = TreeProperties<InnerExpr>::has_differentiable_sub_tree  };
    };




    template<class Tensor, int Order>
    class TreeProperties< Wrapper<Tensor, Order> > {
    public:
        enum { greatest_tensor_order = Order   };
        enum { smallest_tensor_order = Order   };
        enum { has_mat_mat_mul 		 = 0	   };
        enum { has_differentiable_sub_tree = 0 };
    };


    template<class Expr>
    constexpr int is_scalar_expression_tree()
    {
        return TreeProperties<Expr>::greatest_tensor_order == 0;
    }

    template<class Derived>
    inline constexpr bool has_mat_mat_mul(const Expression<Derived> &)
    {
        return TreeProperties<Derived>::has_mat_mat_mul;
    }

}

#endif //UTOPIA_TREE_PROPERTIES_HPP

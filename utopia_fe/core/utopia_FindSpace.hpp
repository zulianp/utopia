#ifndef UTOPIA_FIND_SPACE_HPP
#define UTOPIA_FIND_SPACE_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_Traverse.hpp"

#include <cassert>
#include <memory>

namespace utopia {

    template <class Space, class Expr>
    class FindTestSpace {
    public:
        template <class Any>
        inline constexpr static int visit(const Any &) {
            return TRAVERSE_CONTINUE;
        }

        template <class T>
        inline int visit(const TestFunction<T> &expr) {
            space_ = expr.space_ptr();
            return TRAVERSE_STOP;
        }

        template <class T, int Backend>
        inline int visit(const TestFunction<FiniteElement<T, Backend>> &expr) {
            space_ = utopia::make_ref(expr.space_ptr()->space());
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const TestFunction<ProductFunctionSpace<T>> &expr) {
            prod_space_ = expr.space_ptr();
            space_ = prod_space_->subspace_ptr(0);
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const LinearIntegrator<T> &expr) {
            space_ = expr.test_ptr();
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const LinearIntegrator<ProductFunctionSpace<T>> &expr) {
            prod_space_ = expr.test_ptr();
            space_ = prod_space_->subspace_ptr(0);
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const BilinearIntegrator<T> &expr) {
            space_ = expr.test_ptr();
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const BilinearIntegrator<ProductFunctionSpace<T>> &expr) {
            prod_space_ = expr.test_ptr();
            space_ = prod_space_->subspace_ptr(0);
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const EquationIntegrator<T> &expr) {
            space_ = expr.test_ptr();
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const EquationIntegrator<ProductFunctionSpace<T>> &expr) {
            prod_space_ = expr.test_ptr();
            space_ = prod_space_->subspace_ptr(0);
            return TRAVERSE_STOP;
        }

        FindTestSpace() : space_(nullptr) {}

        inline bool found() const { return static_cast<bool>(space_); }

        template <class ExprTree>
        inline std::shared_ptr<Space> apply(const ExprTree &expr) {
            space_ = nullptr;
            traverse(expr, *this);
            return space_;
        }

        std::shared_ptr<Space> space() const { return space_; }

        std::shared_ptr<ProductFunctionSpace<Space>> prod_space() const { return prod_space_; }

        std::shared_ptr<Space> space_;
        std::shared_ptr<ProductFunctionSpace<Space>> prod_space_;
    };

    template <class Space, class Expr>
    class FindTestSpace<ProductFunctionSpace<Space>, Expr> {
    public:
        FindTestSpace<Space, Expr> find_;

        std::shared_ptr<ProductFunctionSpace<Space>> space() const { return find_.prod_space(); }
    };

    template <class Space, class Expr>
    class FindTrialSpace {
    public:
        template <class Any>
        inline constexpr static int visit(const Any &) {
            return TRAVERSE_CONTINUE;
        }

        template <class T>
        inline int visit(const TrialFunction<T> &expr) {
            space_ = expr.space_ptr();
            return TRAVERSE_STOP;
        }

        template <class T, int Backend>
        inline int visit(const TrialFunction<FiniteElement<T, Backend>> &expr) {
            space_ = utopia::make_ref(expr.space_ptr()->space());
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const TrialFunction<ProductFunctionSpace<T>> &expr) {
            prod_space_ = expr.space_ptr();
            space_ = prod_space_->subspace_ptr(0);
            return TRAVERSE_STOP;
        }

        template <class Coefficient, class Function>
        inline int visit(const Interpolate<Coefficient, Function> &expr) {
            return visit(expr.fun());
        }

        template <class T>
        inline int visit(const BilinearIntegrator<T> &expr) {
            space_ = expr.trial_ptr();
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const BilinearIntegrator<ProductFunctionSpace<T>> &expr) {
            prod_space_ = expr.trial_ptr();
            space_ = prod_space_->subspace_ptr(0);
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const EquationIntegrator<T> &expr) {
            space_ = expr.trial_ptr();
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const EquationIntegrator<ProductFunctionSpace<T>> &expr) {
            prod_space_ = expr.trial_ptr();
            space_ = prod_space_->subspace_ptr(0);
            return TRAVERSE_STOP;
        }

        FindTrialSpace() : space_(nullptr) {}

        inline bool found() const { return static_cast<bool>(space_); }

        template <class ExprTree>
        inline std::shared_ptr<Space> apply(const ExprTree &expr) {
            space_ = nullptr;
            traverse(expr, *this);
            return space_;
        }

        std::shared_ptr<Space> space() const { return space_; }

        std::shared_ptr<ProductFunctionSpace<Space>> prod_space() const { return prod_space_; }

        std::shared_ptr<Space> space_;
        std::shared_ptr<ProductFunctionSpace<Space>> prod_space_;
    };

    template <class Space, class Expr>
    class FindTrialSpace<ProductFunctionSpace<Space>, Expr> {
    public:
        FindTrialSpace<Space, Expr> find_;

        std::shared_ptr<ProductFunctionSpace<Space>> space() const { return find_.prod_space(); }
    };

    template <class Mesh, class Expr>
    class FindMesh {
    public:
        template <class Any>
        inline constexpr static int visit(const Any &) {
            return TRAVERSE_CONTINUE;
        }

        template <class T>
        inline int visit(const TestFunction<T> &expr) {
            mesh_ = &expr.space_ptr()->mesh();
            return TRAVERSE_STOP;
        }

        template <class T>
        inline int visit(const TestFunction<ProductFunctionSpace<T>> &expr) {
            mesh_ = &expr.space_ptr()->subspace(0).mesh();
            return TRAVERSE_STOP;
        }

        FindMesh() : mesh_(nullptr) {}

        inline bool found() const { return mesh_; }

        template <class ExprTree>
        inline Mesh *apply(const ExprTree &expr) {
            mesh_ = nullptr;
            traverse(expr, *this);
            return mesh_;
        }

        Mesh *mesh_;
    };

    template <class Space, class Expr>
    inline auto find_test_space(const Expr &tree) -> std::shared_ptr<Space> {
        FindTestSpace<Space, Expr> fm;
        auto space_ptr = fm.apply(tree);
        return space_ptr;
    }

    template <class Space, class Expr>
    inline auto find_trial_space(const Expr &tree) -> std::shared_ptr<Space> {
        FindTrialSpace<Space, Expr> fm;
        auto space_ptr = fm.apply(tree);
        assert(space_ptr);
        return space_ptr;
    }

    template <class Space, class Expr>
    inline auto find_space(const Expr &tree) -> Space & {
        {
            FindTestSpace<Space, Expr> fm;
            auto space_ptr = fm.apply(tree);

            if (space_ptr) {
                return *space_ptr;
            }
        }

        {
            FindTrialSpace<Space, Expr> fm;
            auto space_ptr = fm.apply(tree);
            assert(space_ptr);
            return *space_ptr;
        }
    }

    template <class Mesh, class Expr>
    inline auto find_mesh(const Expr &tree) -> const Mesh & {
        FindMesh<Mesh, Expr> fm;
        auto mesh_ptr = fm.apply(tree);
        assert(mesh_ptr);
        return *mesh_ptr;
    }

}  // namespace utopia

#endif  // UTOPIA_FIND_SPACE_HPP

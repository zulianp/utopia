#ifndef UTOPIA_LIBMESH_ASSEMBLY_VALUES_HPP
#define UTOPIA_LIBMESH_ASSEMBLY_VALUES_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_Types.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_FunctionalTraits.hpp"
#include "utopia_libmesh_Utils.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_Traverse.hpp"
#include "utopia_SymbolicFunction.hpp"
#include "utopia_ProductFunctionSpace.hpp"
#include "utopia_libmesh_TreeNavigator.hpp"
#include "utopia_DualBasis.hpp"

#include "utopia_fe_core.hpp"

#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/fe_interface.h"

#include "utopia_libmesh_VectorElement.hpp"

namespace utopia {

    template<class T> class Normal;

    class LibMeshAssemblyValues {
    public:
        typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
        typedef TraitsT::FE FE;
        typedef TraitsT::Matrix Matrix;
        typedef TraitsT::Vector Vector;
        typedef TraitsT::DXType DXType;
        typedef TraitsT::JacobianType JacobianType;

        inline std::vector< std::unique_ptr<FE> > &fe()
        {
            return fe_;
        }

        inline const std::vector< std::unique_ptr<FE> > &fe() const
        {
            return fe_;
        }

        inline const std::vector< std::shared_ptr<VectorElement> > &vector_fe() const
        {
            return vector_fe_;
        }

        inline std::vector< std::shared_ptr<VectorElement> > &vector_fe()
        {
            return vector_fe_;
        }

        inline std::vector< std::unique_ptr<FE> > &test()
        {
            return fe_;
        }

        inline std::vector< std::unique_ptr<FE> > &trial()
        {
            return fe_;
        }

        inline const std::vector< std::unique_ptr<FE> > &test()	const
        {
            return fe_;
        }

        // inline const libMesh::DenseMatrix<double> &biorth_weights(const int subspace_id) const
        // {
        //     if(biorth_weights_.size() == 1) {
        //         return biorth_weights_[0];
        //     }

        //     return biorth_weights_[subspace_id];
        // }

        inline const std::vector< std::unique_ptr<FE> > &trial() const
        {
            return fe_;
        }

        inline std::shared_ptr<libMesh::QBase> quad_test() const
        {
            return quad_test_;
        }

        inline std::shared_ptr<libMesh::QBase> quad_trial() const
        {
            if(!quad_trial_) return quad_test_;
            return quad_trial_;
        }

        template<class Expr>
        void init_fe_from(const Expr &expr, const std::shared_ptr<libMesh::QBase> &quad)
        {
            auto space_ptr = find_any_space(expr); assert(bool(space_ptr));
            space_ptr->initialize();
            quadrature_order_ = functional_order(expr, *this);
            const int dim = space_ptr->mesh().mesh_dimension();
            const libMesh::Elem * elem = space_ptr->mesh().elem(current_element_);

            if(is_quad(elem->type())) {
                const int temp = quadrature_order_/2;
                quadrature_order_ = (temp + 1) * 2;
            } else if(is_hex(elem->type())) {
                const int temp = quadrature_order_/2;
                quadrature_order_ = (temp + 2) * 2;
            }

            if(quad) {
                set_up_quadrature(quad);
            } else {
                set_up_quadrature(dim, quadrature_order_);
            }
            block_id_ = elem->subdomain_id();

            const auto &eq_sys = space_ptr->equation_system();
            const std::size_t n_vars = eq_sys.n_vars();
            fe_.resize(n_vars);

            for(std::size_t i = 0; i < n_vars; ++i) {
                auto fe = libMesh::FEBase::build(dim, eq_sys.get_dof_map().variable_type(i));
                fe->attach_quadrature_rule(quad_test().get());
                fe_[i] = std::move(fe);

            }

            //FIXME find-out if this is needed
            dual_fe_.resize(n_vars);

            init_fe_flags(expr);

            for(std::size_t i = 0; i < n_vars; ++i) {
                fe_[i]->reinit(elem);
            }

            //vfe
            for(auto &v_fe_ptr : vector_fe_) {
                if(v_fe_ptr) {
                    v_fe_ptr->init(fe_);
                }
            }
        }

        template<class Expr>
        void reinit_fe_from(const Expr &expr, const std::shared_ptr<libMesh::QBase> &quad)
        {
            auto space_ptr = find_any_space(expr);
            const libMesh::Elem * elem = space_ptr->mesh().elem(current_element_);
            block_id_ = elem->subdomain_id();

            if(quad) {
                set_up_quadrature(quad);
            }

            on_boundary_ = elem->on_boundary();

            const auto &eq_sys = space_ptr->equation_system();
            const std::size_t n_vars = eq_sys.n_vars();

            for(std::size_t i = 0; i < n_vars; ++i) {
                fe_[i]->reinit(elem);
            }

            //vfe
            for(auto &v_fe_ptr : vector_fe_) {
                if(v_fe_ptr) {
                    v_fe_ptr->reinit(fe_);
                }
            }
        }

        inline bool on_boundary() const
        {
            return on_boundary_;
        }

        inline static bool subspaces_are_equal(const libMesh::DofMap &dof_map)
        {
            //FIXME check all info
            const std::size_t n_vars = dof_map.n_variables();
            int order = dof_map.variable_type(0).order;
            for(std::size_t i = 0; i < n_vars; ++i) {
                if(order != dof_map.variable_type(i).order) {
                    return false;
                }
            }

            return true;
        }

        void init_dual(const libMesh::ElemType type)
        {
            if(!has_dual_) {
                return;
            }

            const auto n_fun = dual_fe_.size();
            assert(n_fun == fe_.size());

            for(uint i = 0; i < n_fun; ++i) {
                if(dual_fe_[i].empty()) {
                    dual_fe_[i].init(type);
                }

                dual_fe_[i].compute_values(*fe_[i]);
            }

            if(!vector_fe_.empty()) {
                for(auto &vf : vector_fe_) {
                    //FIXME
                    vf->make_dual(dual_fe_[0].weights_);
                }
            }

            //IMPLEMENT ME
            assert(false);
        }

        template<class Expr>
        void init_side_fe_from(const Expr &expr, const int side)
        {
            auto space_ptr = find_any_space(expr);
            space_ptr->initialize();
            quadrature_order_ = functional_order(expr, *this);
            const int dim = space_ptr->mesh().mesh_dimension();
            const libMesh::Elem * elem = space_ptr->mesh().elem(current_element_);

            auto s_type = side_type(elem->type());

            if(is_quad(s_type)) {
                const int temp = quadrature_order_/2;
                quadrature_order_ = (temp + 1) * 2;
            } else if(is_hex(s_type)) {
                const int temp = quadrature_order_/2;
                quadrature_order_ = (temp + 2) * 2;
            }

            set_up_quadrature(dim-1, quadrature_order_);
            // block_id_ = elem->side_ptr(side)->subdomain_id();
            block_id_ = space_ptr->mesh().get_boundary_info().boundary_id(elem, side);

            const auto &eq_sys = space_ptr->equation_system();
            const std::size_t n_vars = eq_sys.n_vars();
            fe_.resize(n_vars);

            for(std::size_t i = 0; i < n_vars; ++i) {
                auto fe = libMesh::FEBase::build(dim, eq_sys.get_dof_map().variable_type(i));
                fe->attach_quadrature_rule(quad_test().get());
                fe_[i] = std::move(fe);
            }

            init_fe_flags(expr);

            for(std::size_t i = 0; i < n_vars; ++i) {
                fe_[i]->reinit(elem, side);
            }

            //vfe
            for(auto &v_fe_ptr : vector_fe_) {
                if(v_fe_ptr) {
                    v_fe_ptr->init(fe_);
                }
            }

            init_dual(s_type);
        }

        template<class Expr>
        void reinit_side_fe_from(const Expr &expr, const int side)
        {
            auto space_ptr = find_any_space(expr);
            const int dim = space_ptr->mesh().mesh_dimension();
            const libMesh::Elem * elem = space_ptr->mesh().elem(current_element_);

            block_id_ = space_ptr->mesh().get_boundary_info().boundary_id(elem, side);

            const auto &eq_sys = space_ptr->equation_system();
            const std::size_t n_vars = eq_sys.n_vars();

            for(std::size_t i = 0; i < n_vars; ++i) {
                fe_[i]->reinit(elem, side);
            }

            //vfe
            for(auto &v_fe_ptr : vector_fe_) {
                if(v_fe_ptr) {
                    v_fe_ptr->reinit(fe_);
                }
            }

            //init_dual(s_type);
        }

        const DXType &dx() const
        {
            return test()[0]->get_JxW();
        }

        inline int block_id() const
        {
            return block_id_;
        }

        inline void set_current_element(const long current_element)
        {
            current_element_ = current_element;
        }

        inline long current_element() const
        {
            return current_element_;
        }


        inline std::size_t n_shape_functions() const
        {
            std::size_t ret = 0;
            for(auto &f_ptr : fe()) {
                if(f_ptr) ret += f_ptr->n_shape_functions();
            }

            return ret;
        }

        void set_up_quadrature(const std::shared_ptr<libMesh::QBase> &quad)
        {
            quad_test_ = quad;
            quad_trial_ = quad;
            reset_quadrature_ = false;
        }

        void set_up_quadrature(const int dim, const int quadrature_order)
        {
            if(reset_quadrature_) {
                if(quadrature_order == 1) {
                    quad_test_  = std::make_shared<libMesh::QGauss>(dim, libMesh::Order(2));
                } else {
                    quad_test_  = std::make_shared<libMesh::QGauss>(dim, libMesh::Order(quadrature_order));
                }

                quad_trial_ = quad_test_;
                reset_quadrature_ = false;
            } else {
                assert((static_cast<bool>(quad_trial_)));
                assert((static_cast<bool>(quad_test_)));
            }
        }

        void set_has_dual(const bool val) {
            has_dual_ = val; 
        }

        LibMeshAssemblyValues()
        : current_element_(0),
          quadrature_order_(2),
           block_id_(0),
           reset_quadrature_(true),
           on_boundary_(false),
           has_dual_(false)
        {}


        std::vector<DualBasis> &dual_fe() { return dual_fe_; }
        const std::vector<DualBasis> &dual_fe()  const { return dual_fe_; }

    private:
        long current_element_;
        long quadrature_order_;
        int block_id_;
        bool reset_quadrature_;
        bool on_boundary_;
        bool has_dual_;


        std::shared_ptr<libMesh::QBase> quad_trial_;
        std::shared_ptr<libMesh::QBase> quad_test_;
        std::vector< std::unique_ptr<FE> > fe_;
        std::vector<DualBasis> dual_fe_;

        //additional precomputed values
        std::vector<std::shared_ptr<VectorElement>> vector_fe_;


        class FEInitializer {
        public:

            FEInitializer(LibMeshAssemblyValues &ctx)
            : ctx(ctx)
            {}

            template<class Any>
            inline constexpr static int visit(const Any &) { return TRAVERSE_CONTINUE; }


            void init_phi(const LibMeshFunctionSpace &s)
            {
                ctx.fe()[s.subspace_id()]->get_phi();
            }

            void init_phi(const ProductFunctionSpace<LibMeshFunctionSpace> &s)
            {
                s.each([&](const int, const LibMeshFunctionSpace &space) {
                    ctx.fe()[space.subspace_id()]->get_phi();
                });
            }


            void init_JxW(const LibMeshFunctionSpace &s)
            {
                if(s.subspace_id() == 0) {
                    //FIXME
                    ctx.fe()[s.subspace_id()]->get_JxW();
                }
            }

            void init_JxW(const ProductFunctionSpace<LibMeshFunctionSpace> &s)
            {
                s.each([&](const int, const LibMeshFunctionSpace &space) {
                    ctx.fe()[space.subspace_id()]->get_JxW();
                });
            }

            void init_xyz()
            {
                //FIXME
                ctx.fe()[0]->get_xyz();
            }

            template<class T>
            inline int visit(const TrialFunction<T> &expr)
            {
                init_phi(*expr.space_ptr());
                return TRAVERSE_CONTINUE;
            }

            template<class T>
            inline int visit(const TestFunction<T> &expr)
            {
                init_JxW(*expr.space_ptr());
                init_phi(*expr.space_ptr());
                return TRAVERSE_CONTINUE;
            }

            template<class T>
            inline int visit(const Dual<TestFunction<T>> &expr)
            {
                init_phi_dual(*expr.expr().space_ptr());
                return TRAVERSE_CONTINUE;
            }

            template<class T>
            inline int visit(const Gradient<Dual<TestFunction<T>>> &expr)
            {
                init_dphi_dual(*expr.expr().space_ptr());
                return TRAVERSE_CONTINUE;
            }

            template<class T>
            inline int visit(const Divergence<Dual<TestFunction<T>>> &expr)
            {
                init_div_dual(*expr.expr().space_ptr());
                return TRAVERSE_CONTINUE;
            }

            void init_dphi(const LibMeshFunctionSpace &s)
            {
                ctx.fe()[s.subspace_id()]->get_dphi();
            }


            void init_normals()
            {
                ctx.fe()[0]->get_normals();
            }

            void init_dphi(const ProductFunctionSpace<LibMeshFunctionSpace> &s)
            {
                s.each([&](const int, const LibMeshFunctionSpace &space) {
                    ctx.fe()[space.subspace_id()]->get_dphi();
                });

                //vfe
                const std::size_t s_id = s.subspace(0).subspace_id();
                if(ctx.vector_fe().size() <= s_id) {
                    ctx.vector_fe().resize(s_id + 1);
                }

                if(!ctx.vector_fe()[s_id]) {
                    ctx.vector_fe()[s_id] = std::make_shared<VectorElement>();
                    ctx.vector_fe()[s_id]->dim = s.subspace(0).mesh().spatial_dimension();
                    ctx.vector_fe()[s_id]->n_vars = s.n_subspaces();
                    ctx.vector_fe()[s_id]->start_var = s_id;
                }

                ctx.vector_fe()[s_id]->grad_flag = true;
            }

            void init_phi_dual(const LibMeshFunctionSpace &s)
            {
                ctx.set_has_dual(true);
                ctx.dual_fe()[s.subspace_id()].compute_phi = true;
                ctx.dual_fe()[s.subspace_id()].order = s.order();
            }

            void init_dphi_dual(const LibMeshFunctionSpace &s)
            {
                ctx.set_has_dual(true);
                ctx.dual_fe()[s.subspace_id()].compute_dphi = true;
                ctx.dual_fe()[s.subspace_id()].order = s.order();
            }

            void init_div_dual(const LibMeshFunctionSpace &s)
            {
                ctx.set_has_dual(true);
                ctx.dual_fe()[s.subspace_id()].compute_divphi = true;
                ctx.dual_fe()[s.subspace_id()].order = s.order();
            }

            void init_phi_dual(const ProductFunctionSpace<LibMeshFunctionSpace> &space)
            {
                ctx.set_has_dual(true);
                space.each([&](const int, const LibMeshFunctionSpace &s) {
                    ctx.dual_fe()[s.subspace_id()].compute_phi = true;
                    ctx.dual_fe()[s.subspace_id()].order = s.order();
                });
            }

            void init_dphi_dual(const ProductFunctionSpace<LibMeshFunctionSpace> &space)
            {
                ctx.set_has_dual(true);
                space.each([&](const int, const LibMeshFunctionSpace &s) {
                    ctx.dual_fe()[s.subspace_id()].compute_dphi = true;
                    ctx.dual_fe()[s.subspace_id()].order = s.order();
                });
            }

            void init_div_dual(const ProductFunctionSpace<LibMeshFunctionSpace> &space)
            {
               ctx.set_has_dual(true);
               space.each([&](const int, const LibMeshFunctionSpace &s) {
                   ctx.dual_fe()[s.subspace_id()].compute_divphi = true;
                   ctx.dual_fe()[s.subspace_id()].order = s.order();
               });
            }

            template<class Expr>
            inline int visit(const Integral<Expr> &expr)
            {
                ctx.fe()[0]->get_JxW();
                return TRAVERSE_CONTINUE;
            }

            //Gradient
            template<template<class> class Function>
            inline int visit(const Gradient<Function<ProductFunctionSpace<LibMeshFunctionSpace>>> &expr)
            {
                init_dphi(*expr.expr().space_ptr());
                return TRAVERSE_CONTINUE;
            }

            template<template<class> class Function>
            inline int visit(const Gradient<Function<LibMeshFunctionSpace>> &expr)
            {
                init_dphi(*expr.expr().space_ptr());
                return TRAVERSE_CONTINUE;
            }

            template<class C, template<class> class Function>
            inline int visit(const Gradient<Interpolate<C, Function<LibMeshFunctionSpace>> > &expr)
            {
                init_dphi(*expr.expr().fun().space_ptr());
                return TRAVERSE_SKIP_SUBTREE;
            }

            template<class C, template<class> class Function>
            inline int visit(const Gradient<Interpolate<C, Function<ProductFunctionSpace<LibMeshFunctionSpace>>> > &expr)
            {
                init_dphi(*expr.expr().fun().space_ptr());
                return TRAVERSE_SKIP_SUBTREE;
            }

            template<class C, template<class> class Function>
            inline int visit(const Interpolate<C, Function<LibMeshFunctionSpace>> &expr)
            {
                init_phi(*expr.fun().space_ptr());
                return TRAVERSE_CONTINUE;
            }

            template<class C, template<class> class Function>
            inline int visit(const Interpolate<C, Function<ProductFunctionSpace<LibMeshFunctionSpace>>> &expr)
            {
                init_phi(*expr.fun().space_ptr());
                return TRAVERSE_CONTINUE;
            }

            //Divergence
            template<template<class> class Function>
            inline int visit(const Divergence<Function<ProductFunctionSpace<LibMeshFunctionSpace>>> &expr)
            {
                init_dphi(*expr.expr().space_ptr());
                return TRAVERSE_CONTINUE;
            }

            template<template<class> class Function>
            inline int visit(const Divergence<Function<LibMeshFunctionSpace>> &expr)
            {
                init_dphi(*expr.expr().space_ptr());
                return TRAVERSE_CONTINUE;
            }

            //Curl
            template<template<class> class Function>
            inline int visit(const Curl<Function<ProductFunctionSpace<LibMeshFunctionSpace>>> &expr)
            {
                init_dphi(*expr.expr().space_ptr());
                return TRAVERSE_CONTINUE;
            }

            template<template<class> class Function>
            inline int visit(const Curl<Function<LibMeshFunctionSpace>> &expr)
            {
                init_dphi(*expr.expr().space_ptr());
                return TRAVERSE_CONTINUE;
            }

            template<class Out, class F>
            inline int visit(const ContextFunction<Out, F> &expr)
            {
                init_xyz();
                return TRAVERSE_CONTINUE;
            }

            template<class Out, class F>
            inline int visit(const ContextFunction<Out, Normal<F>> &expr)
            {
                init_normals();
                return TRAVERSE_CONTINUE;
            }

#ifdef WITH_TINY_EXPR
            inline int visit(const SymbolicFunction &expr)
            {
                init_xyz();
                return TRAVERSE_CONTINUE;
            }
#endif //WITH_TINY_EXPR

            template<class Expr>
            void apply(const Expr &expr)
            {
                traverse(expr, *this);
            }

            LibMeshAssemblyValues &ctx;
        };

        template<class Expr>
        void init_fe_flags(const Expr &expr)
        {
            FEInitializer fe_init(*this);
            fe_init.apply(expr);
        }

    };

}

#endif //UTOPIA_LIBMESH_ASSEMBLY_VALUES_HPP

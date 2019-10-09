#ifndef UTOPIA_FLOW_WITH_FRACTURES_HPP
#define UTOPIA_FLOW_WITH_FRACTURES_HPP

#include "utopia_Model.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_Integral.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_Flow.hpp"
#include "utopia_LibMeshToMoonolithConvertions.hpp"
#include "moonolith_l2_assembler.hpp"

#include <iostream>

namespace utopia {

    template<class FunctionSpace, class Matrix, class Vector>
    class FlowWithFractures final : public Model<Matrix, Vector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);
        typedef utopia::Traits<FunctionSpace> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;
        typedef typename TraitsT::Vector ElementVector;

        FlowWithFractures(FunctionSpace &space)
        : flow_(space), space_(space), rescale_(1.0)
        {}

        class Fracture : public Configurable {
        public:
            Scalar angle;
            Scalar aperture;
            Scalar length;
            Scalar permeability;

            ElementVector center;
            ElementVector normal;
            ElementMatrix tangents;

            Fracture(const SizeType dim) 
            : angle(1.0), aperture(1.0), length(1.0), permeability(1.0)
            {
                center   = zeros(dim);
                normal   = zeros(dim);
                tangents = zeros(dim-1, dim);
            }

            static std::string coord_str(const SizeType d)
            {
                static const std::vector<std::string> str = {
                   "x", "y", "z", "t"
                };

                return str[d];
            }

            void read(Input &in) override
            {
                in.get("angle", angle);
                in.get("aperture", aperture);
                in.get("length", length);
                in.get("permeability", permeability);

                const SizeType n = size(center); 

                for(SizeType i = 0; i < n; ++i) {
                    Scalar v = 0.0;
                    in.get("center-" + coord_str(i), v);
                    center.set(i, v);

                    // v = 0.0;
                    // in.get("normal-" + coord_str(i), v);
                    // normal.set(i, v);

                    for(SizeType k = 0; k < n-1; ++k) {
                        v = 0.0;
                        in.get("tangent-" + std::to_string(k) + "-" + coord_str(i), v);
                        tangents.set(k, i, v);
                    }
                }

                if(n == 2) {
                    normal.set(0, -tangents.get(0, 1));
                    normal.set(1,  tangents.get(0, 0));
                } else 
                // if(n == 3) 
                {
                    //TODO
                    assert(false);
                    // ElementVector u = zeros(2), v = zeros(2);

                    // u.set(0, tangents.get(0, 1));
                    // u.set(1, tangents.get(0, 0));

                    // v.set(0, tangents.get(1, 1));
                    // v.set(1, tangents.get(1, 0));

                    // normal = cross(u, v);
                    // normal /= norm2(normal);
                }
            }
        };

        class FractureSampler final : public UIFunction<double> {
        public:
        
            Scalar eval(const std::vector<Scalar> &) const override
            {
                const auto &f = *fractures[active_fracture];
                return f.aperture * f.permeability;
            }

            FractureSampler() : active_fracture(0) {}

            std::size_t active_fracture;
            moonolith::Storage<std::shared_ptr<Fracture>> fractures;

        };

        template<int Dim>
        class Cut {};

        template<>
        class Cut<2> {
        public:
            static const int Dim = 2;

            using CellElem = moonolith::Elem<double, Dim, Dim>;
            using SurfElem = moonolith::Edge1<double, Dim>;
            using QuadratureT = moonolith::Quadrature<double, Dim>;
            using Algo = moonolith::BuildQuadratureAlgo<double, Dim, Dim, Dim -1>;

            moonolith::Storage<std::shared_ptr<Fracture>> fractures;
            

            std::shared_ptr<CellElem> cell_elem;
            std::shared_ptr<SurfElem> surf_elem;

            std::shared_ptr<moonolith::Transform<double, Dim, Dim>>         cell_trafo;
            std::shared_ptr<moonolith::AffineTransform<double, Dim-1, Dim>> surf_trafo;

            Algo algo;

            Cut()
            {
               surf_trafo = std::make_shared<moonolith::AffineTransform<double, Dim-1, Dim>>();
            }

            bool compute(
                const int order,
                const libMesh::Elem &elem,
                const libMesh::FEType &type,
                const moonolith::Storage<std::shared_ptr<Fracture>> &fractures,
                moonolith::Storage<std::shared_ptr<libMesh::QBase>> &quadrature
            )
            {
                moonolith::Gauss::get(order, algo.q_rule);

                make(elem, type, cell_elem);
                std::size_t n = fractures.size();
                quadrature.resize(n);

                if(!surf_elem) {
                    surf_elem = std::make_shared<SurfElem>();
                }

                bool intersected = false;
                for(std::size_t i = 0; i < n; ++i) {
                    make_elem(*fractures[i], *surf_elem);

                    if(intersect(*cell_elem, *surf_elem, quadrature[i])) {
                        intersected = true;
                    }

                }
    
                return intersected;
            }

            void make_elem(const Fracture &f, SurfElem &e) const
            {
                e.point(0)[0] = f.center.get(0) - f.tangents.get(0, 0) * (f.length/2.);
                e.point(0)[1] = f.center.get(1) - f.tangents.get(0, 1) * (f.length/2.);
                e.point(1)[0] = f.center.get(0) + f.tangents.get(0, 0) * (f.length/2.);
                e.point(1)[1] = f.center.get(1) + f.tangents.get(0, 1) * (f.length/2.);

                e.node(0)[0] =  f.center.get(0);
                e.node(0)[1] =  f.center.get(1);
            }

            bool intersect(const CellElem &cell, const SurfElem &surf, std::shared_ptr<libMesh::QBase> &q)
            {
                make(cell, algo.master);
                make(surf, algo.slave);

                make_transform(cell, cell_trafo);
                make_transform(surf, *surf_trafo);

                algo.trafo_master = cell_trafo;
                algo.trafo_slave  = surf_trafo;

                if(algo.compute()) {
                    auto q_mortar = std::make_shared<QMortar>(Dim);

                    utopia::convert(
                        algo.q_master,
                        1.0,
                        *q_mortar);

                    q = q_mortar;
                    return true;
                } else {
                    q = nullptr;
                    return false;
                }
            }
        };

        inline bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        {  
            if(!flow_.assemble_hessian_and_gradient(x, hessian, gradient)) {
                return false;
            }

            auto u = trial(space_);
            auto v = test(space_);

            auto bilinear_form = inner(
                    grad(u), 
                    ctx_fun(sampler()) * grad(v)
                ) * dX;

            return assemble(bilinear_form, hessian, true);
        }

        template<class BilinearForm>
        bool assemble(BilinearForm &&bilinear_form, Matrix &hessian, const bool append_mode = true)
        {
            Chrono c;
            c.start();

            const auto &dof_map = space_.dof_map();

            if(!append_mode) {
                LibMeshAssembler::allocate_matrix(dof_map, hessian);
            }

            Cut<2> cut;
            moonolith::Storage<std::shared_ptr<libMesh::QBase>> q;
            USerialMatrix el_mat;
            std::vector<libMesh::dof_id_type> dof_indices;

            AssemblyContext<LIBMESH_TAG> ctx;
            const int order = functional_order(bilinear_form, ctx);
            for(auto e_it = space_.mesh().active_local_elements_begin(); e_it != space_.mesh().active_local_elements_end(); ++e_it) {
                if(!cut.compute(
                    order,
                    **e_it,
                    space_.type(), 
                    fracture_sampler_.fractures,
                    q)) continue;

                ctx.set_current_element((*e_it)->id());

                //reset the fracture id
                fracture_sampler_.active_fracture = 0;
                for(auto q_it = q.begin(); q_it != q.end(); ++q_it)
                {       
                    if(*q_it) {
                        ctx.set_has_assembled(false);   
                        ctx.init( bilinear_form, *q_it );

                        el_mat.set(0.0);

                        FormEvaluator<LIBMESH_TAG> eval;
                        eval.eval(bilinear_form, el_mat, ctx, true);

                        if(ctx.has_assembled()) {
                            dof_map.dof_indices((*e_it), dof_indices);
                            add_matrix(el_mat, dof_indices, dof_indices, hessian);
                        }
                    }

                    //increment the fracture id
                    ++fracture_sampler_.active_fracture;
                }
            }
            
            if(rescale_ != 1.0) {
                hessian *= rescale_;
            }

            c.stop();
            return true;
        }

        inline bool is_linear() const override { return true; }

        inline void clear() override {}

        inline void read(Input &in) override 
        {
            int dim = space_.mesh().spatial_dimension();
            flow_.read(in);
            
            in.get("fractures", [this, dim](Input &in) {
                in.get_all([this, dim](Input &in) {
                    auto f = utopia::make_unique<Fracture>(dim);
                    f->read(in);
                    fracture_sampler_.fractures.push_back(std::move(f));
                });
            });
        }

        inline void rescale(const Scalar rescale)
        {
            rescale_ = rescale;
        }

        inline std::shared_ptr<UIFunction<Scalar>> sampler()
        {
            return make_ref(fracture_sampler_);
        }

    private:
        Flow<FunctionSpace, Matrix, Vector> flow_;
        FunctionSpace &space_;
        Scalar rescale_;
        
        moonolith::Storage<std::shared_ptr<Cut<2>>> intersections_;
        FractureSampler fracture_sampler_;
    };

}

#endif //UTOPIA_FLOW_WITH_FRACTURES_HPP

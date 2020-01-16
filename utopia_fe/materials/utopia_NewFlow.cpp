#include "utopia_NewFlow.hpp"
#include "utopia_LocalMaterial.hpp"

#include "utopia_FEEval_Local.hpp"
#include "utopia_FEEval_MultiTensor.hpp"
#include "utopia_LocalMaterial.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_AssembleLocal.hpp"


namespace utopia {

	template<class Space, typename Scalar>
	void eval(
		const FiniteElement<Space> &fe,
		UIScalarFunction<Scalar> &fun,
		MultiScalard &result)
	{
		const auto &ctx = fe.ctx();
	    fun.set_current_block(ctx.block_id());
	    const auto &pts = ctx.fe()[0]->get_xyz();
	    const auto n = pts.size();

	    result.resize(n);

	    for(std::size_t i = 0; i < n; ++i) {
	        std::vector<Scalar> p = { pts[i](0), pts[i](1), pts[i](2) };
	        result[i] = fun.eval(p);
	    }
	}

	template<class FE>
	class LocalFlow final : public LocalMaterial<FE> {
	public:
		using Scalar   = typename Traits<FE>::Scalar;
		// using SizeType = typename Traits<FE>::SizeType;

		LocalFlow(
			UIScalarFunction<Scalar> &permeability,
			const USerialMatrix &diffusion_tensor)
		: permeability(permeability), diffusion_tensor(diffusion_tensor)
		{}

	    void init(FE &element)
	    {
	    	//Symbolic
	    	auto u = trial(element);
	    	auto subdomain_id = element.ctx().block_id();

	    	 //Symbolic to Numeric
	    	g  = grad(u);
	    	dx = measure(element);

	    	eval(element, permeability, perm);
	    }
	   
	    void assemble(FE &element, USerialMatrix &mat)
	    {
	    	const SizeType n_funs = g.size();

	    	loop(dx.size(), [&](const SizeType &q) {
	    		for(SizeType i = 0; i < n_funs; ++i) {
	    			Ag = diffusion_tensor * g[i][q];
	    			Ag *= (perm[q] * dx[q]);

	    			mat.add(i, i, dot( Ag, g[i][q] ) );

	    			for(SizeType j = i + 1; j < n_funs; ++j) {
	    				const Scalar v = dot( Ag, g[j][q] );
	    				
	    				mat.add(i, j, v);
	    				mat.add(j, i, v);
	    			}
	    		}
	    	});
	    }

	    void assemble(FE &, USerialVector &)
	    {

	    }

	private:
		UIScalarFunction<Scalar> &permeability;
		const USerialMatrix &diffusion_tensor;

		///////////////// BUFFERS //////////////////////
		FormVectord g;
		MultiScalard perm, dx;

		USerialVector Ag;

	};

    template<class FunctionSpace, class Matrix, class Vector>
    NewFlow<FunctionSpace, Matrix, Vector>::NewFlow(FunctionSpace &space)
    : space_(space), forcing_function_(space), rescale_(1.0)
    {}

    template<class FunctionSpace, class Matrix, class Vector>
    bool NewFlow<FunctionSpace, Matrix, Vector>::assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient)
    {
    	using FE = utopia::FiniteElement<FunctionSpace>;

        Chrono c;
        c.start();

        LocalFlow<FE> assembler(permeability_, diffusion_tensor_);

        assemble(
            space_,
            //FIXME
            ctx_fun(permeability_.sampler()) * inner(grad(trial(space_)), grad(trial(space_))) * dX,
            hessian,
            [&](FE &element, USerialMatrix &mat)
            {
                assembler.init(element);
                assembler.assemble(element, mat);
            }
        );


        forcing_function_.eval(x, gradient);

        //for newton methods
        // gradient -= hessian * x;

        if(rescale_ != 1.0) {
            hessian *= rescale_;
            gradient *= rescale_;
        }

        c.stop();
        std::cout << "Flow assemly time: " << c << std::endl;

        assemble_lower_dimensional_features(x, hessian, gradient);

        const Scalar sum_g = sum(gradient);
        std::cout << "sum_g: " << sum_g << std::endl;


        return true;
    }

    template<class FunctionSpace, class Matrix, class Vector>
    bool NewFlow<FunctionSpace, Matrix, Vector>::assemble_lower_dimensional_features(const Vector &x, Matrix &hessian, Vector &gradient)
    {
        if(lower_dimensional_tags_.empty()) {
            return true;
        }

        Chrono c;
        c.start();

        auto u = trial(space_);
        auto v = test(space_);

        const std::size_t n = lower_dimensional_tags_.size();

        Matrix trace_hessian;
        for(std::size_t i = 0; i < n; ++i) {
            auto side = lower_dimensional_tags_[i];

            auto bilinear_form = surface_integral( 
                inner(
                /* diffusion_tensor_ * */ 
                grad(u), 
                ctx_fun(lower_dimensional_permeability_[i]) * grad(v)),
                side
            );

            utopia::assemble(bilinear_form, trace_hessian);

            if(rescale_ != 1.0) {
                trace_hessian *= rescale_;
            }

            hessian += trace_hessian;
        }


        c.stop();
        std::cout << "Flow assemly (sub-dimensional) time: " << c << std::endl;
        return false;
    }

    template<class FunctionSpace, class Matrix, class Vector>
    void NewFlow<FunctionSpace, Matrix, Vector>::read(Input &in) 
    {
        read_permeability_tensor(in);
        in.get("permeability-function", permeability_);
        in.get("forcing-function", forcing_function_);

        lower_dimensional_tags_.clear();
        lower_dimensional_tags_.clear();

        std::cout << "lower-dimensional-permeability:\n";

        in.get("lower-dimensional-permeability", [this](Input &in) {
            in.get_all([this](Input &in) {
                int tag = -1;

                Scalar value = 1.0;
                in.get("value", value);
                in.get("side", tag);

                std::cout << "side(" << tag << "): " << value << std::endl;
                
                if(tag != -1) {
                    auto fun = std::make_shared<UIConstantFunction<Scalar>>(value);

                    lower_dimensional_permeability_.push_back(fun);
                    lower_dimensional_tags_.push_back(tag);
                } else {
                    std::cerr << "[Error] malformed input for lower-dimensional-permeability!" << std::endl;
                }
            });
        });
    }


    template<class FunctionSpace, class Matrix, class Vector>
    void NewFlow<FunctionSpace, Matrix, Vector>::read_permeability_tensor(Input &in)
    {
        Scalar constant_permeability = 1.;
        in.get("permeability",   constant_permeability);

        Scalar permeabilities[3] = {
            constant_permeability,
            constant_permeability,
            constant_permeability
        };
       
        in.get("permeability-x", permeabilities[0]);
        in.get("permeability-y", permeabilities[1]);
        in.get("permeability-z", permeabilities[2]);

        const int dim = space_.mesh().spatial_dimension();

        diffusion_tensor_ = identity(dim, dim);

        {
            Write<ElementMatrix> w(diffusion_tensor_);
            for(int i = 0; i < dim; ++i) {
                diffusion_tensor_.set(i, i, permeabilities[i]);
            }
        }

        std::cout << "global permeabilty tensor: ";

        for(int i = 0; i < dim; ++i) {
            std::cout << permeabilities[i] << " ";
        }

        std::cout << std::endl;
    }

	template class NewFlow<LibMeshFunctionSpace, USparseMatrix, UVector>;
}
#ifndef UTOPIA_UI_SCALAR_SAMPLER_HPP
#define UTOPIA_UI_SCALAR_SAMPLER_HPP

#include "utopia_ui.hpp"
#include "utopia_CSV.hpp"
#include "utopia_FEFunction.hpp"

#include <fstream>
#include <cassert>
#include <cstdio>

#include <libmesh/exact_solution.h>

namespace utopia {

	template<typename Scalar>
	class UIFunction {
	public:
		virtual ~UIFunction() {}
		virtual Scalar eval(const std::vector<Scalar> &x) const = 0;
	};

	template<typename Scalar>
	class UIConstantFunction final : public UIFunction<Scalar> {
	public:
		virtual ~UIConstantFunction() {}

		UIConstantFunction(const Scalar val)
		: val_(val)
		{}

		inline Scalar eval(const std::vector<Scalar> &) const
		{
			return val_;
		}

		const Scalar &value() const
		{
			return val_;
		}

	private:
		Scalar val_;
	};

	template<typename Scalar>
	class UIBoxedFunction final : public UIFunction<Scalar> {
	public:

		UIBoxedFunction(
			const std::vector<Scalar> &lowbo,
			const std::vector<Scalar> &upbo,
			const std::shared_ptr<UIFunction<Scalar>> &fun
			)
		: lowbo_(lowbo), upbo_(upbo), fun_(fun)
		{}

		inline Scalar eval(const std::vector<Scalar> &x) const override
		{
			const auto n = std::min(x.size(), lowbo_.size());
			for(std::size_t i = 0; i < n; ++i) {
				if(lowbo_[i] > x[i] || upbo_[i] < x[i]) {
					return 0.;
				}
			}

			auto value = fun_->eval(x);
			return value;
		}

	private:
		std::vector<Scalar> lowbo_;
		std::vector<Scalar> upbo_;
		std::shared_ptr<UIFunction<Scalar>> fun_;
	};

	template<typename Scalar>
	class UILambdaFunction final : public UIFunction<Scalar> {
	public:
		using F = std::function<Scalar(const std::vector<Scalar> &)>;

		UILambdaFunction(F fun)
		: fun_(fun)
		{}

		inline Scalar eval(const std::vector<Scalar> &x) const override
		{
			return fun_(x);
		}

	private:
		F fun_;
	};

	template<typename F>
	std::shared_ptr<UIFunction<double>> lambda_fun(
		F fun
	) {	
		return std::make_shared<UILambdaFunction<double>>(fun);
	}

	template<typename F>
	std::shared_ptr<UIFunction<double>> boxed_fun(
		const std::vector<double> &lowbo,
		const std::vector<double> &upbo,
		F fun
	) {	
		return std::make_shared<UIBoxedFunction<double>>(lowbo, upbo, lambda_fun(fun));
	}
	

	template<typename Scalar_>
	class ContextFunction<std::vector<Scalar_>, UIFunction<Scalar_>> : public Expression< ContextFunction<std::vector<Scalar_>, UIFunction<Scalar_>> >{
	public:
		static const int Order = 0;
		typedef Scalar_ Scalar;

		ContextFunction(const std::shared_ptr<UIFunction<Scalar>> &fun)
		: fun_(fun)
		{}

		template<int Backend>
		auto eval(const AssemblyContext<Backend> &ctx) const -> std::vector<Scalar>
		{
			const auto &pts = ctx.fe()[0]->get_xyz();

			const auto n = pts.size();
			std::vector<Scalar> ret(n);

			for(std::size_t i = 0; i < n; ++i) {
				std::vector<Scalar> p = { pts[i](0), pts[i](1), pts[i](2) };
				ret[i] = fun_->eval(p);
			}

			return ret;
		}

	private:
		std::shared_ptr<UIFunction<Scalar>> fun_;
	};

	template<typename Scalar>
	inline ContextFunction<std::vector<Scalar>, UIFunction<Scalar> > ctx_fun(const std::shared_ptr<UIFunction<Scalar>> &fun)
	{
		return ContextFunction<std::vector<Scalar>, UIFunction<Scalar> >(fun);
	}

	template<typename Scalar>
	class UIScalarSampler final : public Configurable, public UIFunction<Scalar> {
	public:

		UIScalarSampler() {}

		~UIScalarSampler() {}

		void read(Input &is) override {
			std::string file = "";
			is.get("file", file);

			std::fill(std::begin(min_), std::end(min_), 0.);
			std::fill(std::begin(max_), std::end(max_), 0.);
			std::fill(std::begin(dims_), std::end(dims_), 0);

			is.get("min-x", min_[0]);
			is.get("min-y", min_[1]);
			is.get("min-z", min_[2]);

			is.get("max-x", max_[0]);
			is.get("max-y", max_[1]);
			is.get("max-z", max_[2]);

			is.get("nx", dims_[0]);
			is.get("ny", dims_[1]);
			is.get("nz", dims_[2]);

			n = 0;

			std::size_t n_values = 1;
			for(std::size_t i = 0; i < 3; ++i) {
				range_[i] = max_[i] - min_[i];

				if(dims_[i] != 0) {
					n++;
					n_values *= dims_[i];
				}
			}

			values_.reserve(n_values);
			std::ifstream ifs(file.c_str());
			
			while(ifs.good()) {
				Scalar val;
				ifs >> val;
				values_.push_back(val);
			}

			ifs.close();

			std::size_t n_read_values = values_.size();
			assert(n_read_values == n_values);
			
			if(n_read_values != n_values) {
				std::cout << 
				"number of values is not consistent" 
				"with the dimensions provided expected " 
				<< n_values << " but found " << values_.size() << std::endl;
			}

		}

		Scalar eval(const std::vector<Scalar> &x) const override
		{
			long ind = index(x);
			if(ind < 0 || ind >= values_.size()) {
				assert(false);
				return 0.;
			}
			// std::cout << ind << " = " << (values_[ind]) << std::endl;
			return values_[ind];
		}

		long index(const std::vector<Scalar> &x) const
		{
			long ret = 0;


			std::size_t offset = 1;
			for(std::size_t i = 0; i < n; ++i) {
				long ind = std::min(
					long(dims_[i] - 1), 
					long(floor( (x[i] - min_[i])/range_[i] * dims_[i] ))
				);


				if(ind < 0 || ind >= dims_[i]) {
					//out of range
					return -1;
				}


				ret += offset * ind;
				offset *= dims_[i];
			}

			assert(std::size_t(ret) < values_.size());
			return ret;
		}

		inline bool empty() const
		{
			return values_.empty();
		}

		void describe(std::ostream &os) const
		{
			for(std::size_t i = 0; i < 3; ++i) {
				os << "[" << min_[i] << " " << max_[i] << "] " << range_[i] << " \n";
				os << "dims(" << i << ") " << dims_[i] << "\n";
			}
		}

	private:

		Scalar   min_[3];
		Scalar   max_[3];
		Scalar   range_[3];
		SizeType dims_[3];
		SizeType n;
		std::vector<Scalar> values_;

	};

	template<typename Scalar>
	class Normal final {};


	template<typename Scalar_>
	class ContextFunction<
		std::vector<libMesh::VectorValue<Scalar_>>,
		Normal<Scalar_>
		> final : 
		public Expression< 
				ContextFunction<std::vector<libMesh::VectorValue<Scalar_>>, Normal<Scalar_>>
				>{
	public:
		static const int Order = 1;
		typedef Scalar_ Scalar;

		ContextFunction()
		{}

		template<int Backend>
		auto eval(const AssemblyContext<Backend> &ctx) const -> std::vector<libMesh::VectorValue<Scalar_>>
		{
			const auto &n = ctx.fe()[0]->get_normals();
			auto nn = n.size();

			std::vector<libMesh::VectorValue<Scalar_>> normals(nn);
			for(std::size_t i = 0; i < nn; ++i) {
				for(int d = 0; d < LIBMESH_DIM; ++d) {
					normals[i](d) = n[i](d);
				}
			}

			return normals;
		}

	};

	inline ContextFunction<std::vector<libMesh::VectorValue<double>>, Normal<double>> normal()
	{
		return ContextFunction<std::vector<libMesh::VectorValue<double>>, Normal<double>>();
	}


}

#endif //UTOPIA_UI_SCALAR_SAMPLER_HPP

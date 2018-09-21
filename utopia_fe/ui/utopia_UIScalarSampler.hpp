#ifndef UTOPIA_UI_SCALAR_SAMPLER_HPP
#define UTOPIA_UI_SCALAR_SAMPLER_HPP

#include "utopia_ui.hpp"
#include "utopia_CSV.hpp"
#include "utopia_FEFunction.hpp"

#include <fstream>
#include <cassert>
#include <cstdio>

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

	private:
		Scalar val_;
	};

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
	class UIScalarSampler final : public Serializable, public UIFunction<Scalar> {
	public:

		UIScalarSampler() {}

		~UIScalarSampler() {}

		void read(InputStream &is) override {
			std::string file = "";
			is.read("file", file);

			std::fill(std::begin(min_), std::end(min_), 0.);
			std::fill(std::begin(max_), std::end(max_), 0.);
			std::fill(std::begin(dims_), std::end(dims_), 0);

			is.read("min-x", min_[0]);
			is.read("min-y", min_[1]);
			is.read("min-z", min_[2]);

			is.read("max-x", max_[0]);
			is.read("max-y", max_[1]);
			is.read("max-z", max_[2]);

			for(std::size_t i = 0; i < 3; ++i) {
				range_[i] = max_[i] - min_[i];
			}

			std::ifstream ifs(file.c_str());
			
			if(ifs.good()) {
				std::string line;
				std::getline(ifs, line);
				n = std::sscanf(line.c_str(), "%ld %ld %ld", &dims_[0], &dims_[1], &dims_[2]); assert(n <= 2);

				//for now
				CSV csv;
				csv.read(ifs);

				csv.write(std::cout);

				dims_[0] = csv.n_rows();
				dims_[1] = csv.n_cols();

				values_.resize(csv.n_rows() * csv.n_cols());
				std::size_t index = 0;

				for(auto i = 0; i < csv.n_rows(); ++i) {
					for(auto j = 0; j < csv.n_cols(); ++j) {
						csv.get(i, j, values_[index++]);
					}	
				}

			}

			ifs.close();

		}

		Scalar eval(const std::vector<Scalar> &x) const override
		{
			long ind = index(x);
			if(ind < 0) {
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

}

#endif //UTOPIA_UI_SCALAR_SAMPLER_HPP

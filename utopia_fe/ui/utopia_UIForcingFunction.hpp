#ifndef UTOPIA_UI_FORCING_FUNCTION_HPP
#define UTOPIA_UI_FORCING_FUNCTION_HPP

#include "utopia_ForcedMaterial.hpp"
#include "utopia_ui.hpp"

namespace utopia {
	template<class FunctionSpace, class Vector>
	class UIForcingFunction final : public CompositeForcingFunction<Vector>, public Configurable {
	public:

		UIForcingFunction(FunctionSpace &V) : V_(V) {}

		~UIForcingFunction() {}

		void read(Input &is) override {
			is.read_all([this](Input &is) {

				int block = -1;

				std::string type = "volume";
				is.read("block", block);
				is.read("type", type);

#ifdef WITH_TINY_EXPR
				std::string value;
				is.read("value", value);
				auto f = symbolic(value);
#else
				double value = 0.;
				is.read("value", value);
				auto f = coeff(value);
#endif //WITH_TINY_EXPR

				if(type == "surface") {
					auto v = test(V_);
					auto l_form = surface_integral(inner(f, v), block);

					auto ff = std::make_shared<ConstantForcingFunction<Vector>>();
					ff->init(l_form);
					this->add(ff);
				} else {
					auto v = test(V_);
					auto l_form = integral(inner(f, v), block);

					auto ff = std::make_shared<ConstantForcingFunction<Vector>>();
					ff->init(l_form);
					this->add(ff);
				}
			});
		}

	private:
		FunctionSpace &V_;
	};


	template<class FunctionSpace, class Vector>
	class UIForcingFunction<ProductFunctionSpace<FunctionSpace>, Vector> final : public CompositeForcingFunction<Vector>, public Configurable {
	public:

		UIForcingFunction(ProductFunctionSpace<FunctionSpace> &V) : V_(V) {}

		virtual ~UIForcingFunction() {}

		void read(Input &is) override {
			is.read_all([this](Input &is) {

				int block = -1;
				int coord = 0;

				std::string type = "volume";
				is.read("block", block);
				is.read("coord", coord);
				is.read("type", type);

#ifdef WITH_TINY_EXPR
				std::string value;
				is.read("value", value);
				auto f = symbolic(value);
#else
				double value = 0.;
				is.read("value", value);
				auto f = coeff(value);
#endif //WITH_TINY_EXPR

				if(type == "surface") {
					auto v = test(V_[coord]);
					auto l_form = surface_integral(inner(f, v), block);

					auto ff = std::make_shared<ConstantForcingFunction<Vector>>();
					ff->init(l_form);
					this->add(ff);
				} else {
					auto v = test(V_[coord]);
					auto l_form = integral(inner(f, v), block);

					auto ff = std::make_shared<ConstantForcingFunction<Vector>>();
					ff->init(l_form);
					this->add(ff);
				}
			});
		}

	private:
		ProductFunctionSpace<FunctionSpace> &V_;
	};
}


#endif //UTOPIA_UI_FORCING_FUNCTION_HPP

#ifndef UTOPIA_SDF_RESAMPLE_HPP
#define UTOPIA_SDF_RESAMPLE_HPP

#include "utopia_Field.hpp"
#include "utopia_Traits.hpp"

#include <memory>


namespace utopia {
	template<class FunctionSpace>
	class SDFResample : public Configurable {
	public:
		using Traits = utopia::Traits<FunctionSpace>;
		using Field = utopia::Field<FunctionSpace>;
		using Vector = typename Traits::Vector;
		using Scalar = typename Traits::Scalar;

		SDFResample();
		~SDFResample();
		void read(Input &in) override;

		bool apply(FunctionSpace &space, Vector &out);

	private:
		class Impl;
		std::unique_ptr<Impl> impl_;
	};
}

#endif //UTOPIA_SDF_RESAMPLE_HPP

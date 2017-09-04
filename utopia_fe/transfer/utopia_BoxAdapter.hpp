#ifndef UTOPIA_BOX_ADAPTER_HPP
#define UTOPIA_BOX_ADAPTER_HPP 


#include "Box.hpp"
#include "moonolith_serializable.hpp"
#include "moonolith_describable.hpp"
#include "moonolith_input_stream.hpp"
#include "moonolith_output_stream.hpp"
#include "moonolith_bounding_volume_with_span.hpp"


namespace utopia {

	template<int Dimension>
	using BoxBoxAdapter = moonolith::AABBWithKDOPSpan<Dimension, double>;

	// template<int Dimension>
	// using BoxBoxAdapter = moonolith::AABBWithSpan<Dimension, double>;
}

#endif //UTOPIA_BOX_ADAPTER_HPP


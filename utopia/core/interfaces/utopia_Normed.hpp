#ifndef UTOPIA_NORMED_HPP
#define UTOPIA_NORMED_HPP

namespace utopia {

	template<typename Scalar_>
	class Normed {
	public:
		using Scalar = Scalar_;

		virtual ~Normed() {}
		virtual Scalar norm_infty() const = 0;
		// virtual norm_infty() const = 0;
	};
}

#endif //UTOPIA_NORMED_HPP
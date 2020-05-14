#ifndef UTOPIA_NORMED_HPP
#define UTOPIA_NORMED_HPP

namespace utopia {

	template<typename Scalar_>
	class Normed {
	public:
		using Scalar = Scalar_;

                virtual ~Normed() = default;
                virtual Scalar norm_infty() const = 0;
		virtual Scalar norm1() const = 0;
		virtual Scalar norm2() const = 0;
	};
}

#endif //UTOPIA_NORMED_HPP
#ifndef UTOPIA_I_FUNCTION_SPACE_HPP
#define UTOPIA_I_FUNCTION_SPACE_HPP

#include <memory>

namespace utopia {

	template<typename Integer, typename Scalar>
	class IElem {
	public:
		virtual ~IElem() {}
		virtual std::shared_ptr<IElem> side(const Integer side_num) = 0;
	};

	template<typename Integer, typename Scalar>
	class IFunctionSpace {
	public:
		using Elem = utopia::IElem<Integer, Scalar>;

		virtual ~IFunctionSpace() {}
		virtual void dofs(std::vector<Integer> &indices) = 0;
		virtual void dofs(const Integer var_num, std::vector<Integer> &indices) = 0;
		virtual Integer n_variables() const = 0;
		virtual Integer n_dofs() const = 0;
		virtual Integer n_local_dofs() const = 0;
		virtual Integer n_elements() const = 0;
		virtual std::unique_ptr<Elem> elem(const Integer id) = 0;
		virtual std::unique_ptr<const Elem> elem(const Integer id) const = 0;
	};
}

#endif //UTOPIA_I_FUNCTION_SPACE_HPP

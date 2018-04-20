#ifndef UTOPIA_CLONABLE_HPP
#define UTOPIA_CLONABLE_HPP 

namespace utopia {
	class Clonable {
	public:
		virtual ~Clonable() {}
		virtual Clonable * clone() const = 0;

	};
}

#endif //UTOPIA_CLONABLE_HPP

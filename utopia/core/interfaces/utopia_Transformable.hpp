#ifndef UTOPIA_TRANSFORMABLE_HPP
#define UTOPIA_TRANSFORMABLE_HPP

#include "utopia_Operators.hpp"

namespace utopia {

	template<typename T>
	class Transformable {
	public:
		virtual ~Transformable() {}

		virtual void transform(const Sqrt &) = 0;
		virtual void transform(const Pow2 &) = 0;
		virtual void transform(const Log &)  = 0;
		virtual void transform(const Exp &)  = 0;
		virtual void transform(const Cos &)  = 0;
		virtual void transform(const Sin &)  = 0;
		virtual void transform(const Abs &)  = 0;
		virtual void transform(const Minus &) = 0;

		virtual void transform(const Pow &p)  = 0;
		virtual void transform(const Reciprocal<T> &f) = 0;

		// //UNARY
		// class Sqrt {
		// public:
		//     std::string get_class() const { return "Sqrt"; }

		//     template<typename T>
		//     inline static T apply(const T &x) {
		//         using std::sqrt;
		//         return sqrt(x);
		//     }
		// };

		// class Pow2 {
		// public:
		//     std::string get_class() const { return "Pow2"; }

		//     template<typename T>
		//     inline static T apply(const T &x) {
		//         return x * x;
		//     }
		// };


		// class Pow {
		// public:
		//     std::string get_class() const { return "Pow"; }

		//     template<typename T>
		//     inline T apply(const T &x) const {
		//         return std::pow(x,a_);
		//     }

		//    Pow(const double &a):a_(a){}
		//    double a_;

		// };

		// class Log {
		// public:
		//     std::string get_class() const { return "Log"; }

		//     template<typename T>
		//     inline static T apply(const T &x) {
		//         using std::log;
		//         return log(x);
		//     }
		// };

		// class Exp {
		// public:
		//     std::string get_class() const { return "Exp"; }

		//     template<typename T>
		//     inline static T apply(const T &x) {
		//         using std::exp;
		//         return exp(x);
		//     }
		// };

		// class Cos {
		// public:
		//     std::string get_class() const { return "Cos"; }

		//     template<typename T>
		//     inline static T apply(const T &x) {
		//         using std::cos;
		//         return cos(x);
		//     }
		// };


		// class Sin {
		// public:
		//     std::string get_class() const { return "Sin"; }

		//     template<typename T>
		//     inline static T apply(const T &x) {
		//         using std::sin;
		//         return sin(x);
		//     }
		// };

		// class Transpose {
		// public:
		//     std::string get_class() const { return "Transpose"; }
		// };

		// class Abs {
		// public:
		//     std::string get_class() const { return "Abs"; }

		//     template<typename T>
		//     inline static T apply(const T &x) {
		//         using std::abs;
		//         return abs(x);
		//     }
		// };


		// template<class Expr>
		// std::string GetClass() {
		//     return Expr().get_class();
		// }

		// //Special (FIND A NAME)
		// template<typename  T>
		// class Reciprocal {
		// public:
		//     typedef T Scalar;

		//     std::string get_class() const { return "Reciprocal"; }

		//     template<typename T2>
		//     inline T2 apply(const T2 &x) const {
		//         return numerator_/x;
		//     }

		//     Reciprocal(const T &numerator)
		//             : numerator_(numerator)
		//     {}

		//     inline const T &numerator() const { return numerator_; }
		// private:
		//     const T numerator_;
		// };

	};
}

#endif //UTOPIA_TRANSFORMABLE_HPP

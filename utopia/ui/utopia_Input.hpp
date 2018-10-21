#ifndef UTOPIA_INPUT_STREAM_HPP
#define UTOPIA_INPUT_STREAM_HPP

#include <memory>
#include <string>
#include <vector>
#include <set>
#include <functional>

#include "utopia_Base.hpp"
#include "utopia_Convertible.hpp"

namespace utopia {
	class Path;
	class Input;

	class Configurable {
	public:
		virtual ~Configurable() {}
		virtual void read(Input &is)        = 0;
		// virtual void write(OutputStream &is) const = 0;
	};

	class Input {
	public:
		virtual ~Input() {}

		// virtual bool open(const Path &path) = 0;

		virtual SizeType size() const = 0;
		virtual void get(std::vector<std::shared_ptr<IConvertible>> &values) = 0;

		// template<class T>
		// void get(std::vector<T> &vec) {
		// 	auto n = size();
		// 	if(n == 0) return;

		// 	vec.resize(n);

		// 	array_start();

		// 	for(auto &v : vec) {
		// 		get(v);
		// 		next();
		// 	}

		// 	array_finish();
		// }

		


		// template<class T>
		// void get(const std::string &name, std::vector<T> &vec) {
		// 	get(name, [&vec](Input &sub_is) {
		// 		sub_is.get(vec);
		// 	});
		// }

		// template<class T>
		// void get(std::set<T> &s) {
		// 	auto n = size();

		// 	array_start();

		// 	for(SizeType i = 0; i < n; ++i) {
		// 		T v;
		// 		get(v);
		// 		s.insert(v);
		// 		next();
		// 	}

		// 	array_finish();
		// }	

		virtual void get_all(std::function<void(Input &)> lambda) = 0;

		// virtual void get(bool &val) = 0;
		// virtual void get(double &val) = 0;
		// virtual void get(int &val) = 0;
		// virtual void get(SizeType &val) = 0;
		// virtual void get(std::string &val) = 0;
		// virtual void get(Configurable &val) = 0;
		// virtual void get(std::function<void(Input &)> lambda) = 0;

		virtual void get(const std::string &key, bool &val) = 0;
		virtual void get(const std::string &key, double &val) = 0;
		virtual void get(const std::string &key, int &val) = 0;
		virtual void get(const std::string &key, SizeType &val) = 0;
		virtual void get(const std::string &key, std::string &val) = 0;
		virtual void get(const std::string &key, Configurable &val) = 0;
		virtual void get(const std::string &key, std::function<void(Input &)> lambda) = 0;

		virtual bool good() const = 0;

	protected:
		// virtual void next()   = 0;
		// virtual void array_start()  = 0;
		// virtual void array_finish() = 0;
	};
}

#endif //UTOPIA_INPUT_STREAM_HPP

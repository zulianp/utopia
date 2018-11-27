#ifndef UTOPIA_INPUT_STREAM_HPP
#define UTOPIA_INPUT_STREAM_HPP

#include <memory>
#include <string>
#include <vector>
#include <set>
#include <functional>
#include <ostream>

#include "utopia_Base.hpp"
#include "utopia_Convertible.hpp"
#include "utopia_Path.hpp"

namespace utopia {
	class Path;
	class Input;

	class Configurable {
	public:
		virtual ~Configurable() {}
		virtual void read(Input &is) = 0;
		virtual void print_usage(std::ostream &os = std::cout) const;
		virtual bool import(const Path &path);
		virtual bool import(
			const std::string &key,
			const Path &path);
	};

	class Input  /* : public Clonable */ {
	public:
		Input() {}
		virtual ~Input() {}

		virtual SizeType size() const = 0;
		virtual void get(std::vector<std::shared_ptr<IConvertible>> &values) = 0;
		virtual void get_all(std::function<void(Input &)> lambda) = 0;

		virtual void get(const std::string &key, bool &val) = 0;
		virtual void get(const std::string &key, double &val) = 0;
		virtual void get(const std::string &key, int &val) = 0;
		virtual void get(const std::string &key, long &val) = 0;
		virtual void get(const std::string &key, unsigned long &val) = 0;
		// virtual void get(const std::string &key, SizeType &val) = 0;
		virtual void get(const std::string &key, std::string &val) = 0;
		virtual void get(const std::string &key, Configurable &val) = 0;
		virtual void get(const std::string &key, std::function<void(Input &)> lambda) = 0;
		virtual bool good() const = 0;

	private:
		Input(const Input &other) {}
		Input &operator=(const Input &other) {
			return *this;
		}
	};
}

#endif //UTOPIA_INPUT_STREAM_HPP

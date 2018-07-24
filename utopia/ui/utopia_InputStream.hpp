#ifndef UTOPIA_INPUT_STREAM_HPP
#define UTOPIA_INPUT_STREAM_HPP

#include <memory>
#include <string>

#include "utopia_Base.hpp"

namespace utopia {
	class Path;

	class InputStream {
	public:
		virtual ~InputStream() {}

		virtual bool open(const Path &path) = 0;
		virtual bool object_begin(const std::string &name) = 0;
		virtual bool object_end() = 0;

		virtual void read(double &val) = 0;
		virtual void read(int &val) = 0;
		virtual void read(SizeType &val) = 0;
		virtual void read(std::string &val) = 0;

		virtual void read(const std::string &key, double &val) = 0;
		virtual void read(const std::string &key, int &val) = 0;
		virtual void read(const std::string &key, SizeType &val) = 0;
		virtual void read(const std::string &key, std::string &val) = 0;
		virtual bool good() const = 0;

	};
}

#endif //UTOPIA_INPUT_STREAM_HPP

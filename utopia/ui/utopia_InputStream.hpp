#ifndef UTOPIA_INPUT_STREAM_HPP
#define UTOPIA_INPUT_STREAM_HPP

#include <memory>
#include <string>
#include <vector>
#include <set>
#include <functional>

#include "utopia_Base.hpp"

namespace utopia {
	class Path;
	class InputStream;

	class Serializable {
	public:
		virtual ~Serializable() {}
		virtual void read(InputStream &is)        = 0;
		// virtual void write(OutputStream &is) const = 0;
	};

	class InputStream {
	public:
		virtual ~InputStream() {}

		virtual bool open(const Path &path) = 0;

		virtual SizeType size() const = 0;

		template<class T>
		void read(std::vector<T> &vec) {
			vec.resize(size());

			array_start();

			for(auto &v : vec) {
				read(v);
				next();
			}

			array_finish();
		}

		void read_all(std::function<void(InputStream &)> lambda) {
			auto n = size();

			array_start();

			for(SizeType i = 0; i < n; ++i) {
				lambda(*this);
				next();
			}

			array_finish();
		}

		template<class T>
		void read(const std::string &name, std::vector<T> &vec) {
			read(name, [&vec](InputStream &sub_is) {
				sub_is.read(vec);
			});
		}

		template<class T>
		void read(std::set<T> &s) {
			auto n = size();

			array_start();

			for(SizeType i = 0; i < n; ++i) {
				T v;
				read(v);
				s.insert(v);
				next();
			}

			array_finish();
		}

		virtual void read(bool &val) = 0;
		virtual void read(double &val) = 0;
		virtual void read(int &val) = 0;
		virtual void read(SizeType &val) = 0;
		virtual void read(std::string &val) = 0;
		virtual void read(Serializable &val) = 0;
		virtual void read(std::function<void(InputStream &)> lambda) = 0;

		virtual void read(const std::string &key, bool &val) = 0;
		virtual void read(const std::string &key, double &val) = 0;
		virtual void read(const std::string &key, int &val) = 0;
		virtual void read(const std::string &key, SizeType &val) = 0;
		virtual void read(const std::string &key, std::string &val) = 0;
		virtual void read(const std::string &key, Serializable &val) = 0;
		virtual void read(const std::string &key, std::function<void(InputStream &)> lambda) = 0;

		virtual bool good() const = 0;

	protected:
		virtual void next()   = 0;
		virtual void array_start()  = 0;
		virtual void array_finish() = 0;
	};
}

#endif //UTOPIA_INPUT_STREAM_HPP

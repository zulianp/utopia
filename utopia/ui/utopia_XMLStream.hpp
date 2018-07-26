#ifndef UTOPIA_XML_STREAM_HPP
#define UTOPIA_XML_STREAM_HPP

#include <memory>
#include "utopia_Base.hpp"
#include "utopia_Path.hpp"
#include "utopia_InputStream.hpp"


namespace utopia {

	class XMLInputStream final : public InputStream {
	public:
		XMLInputStream();
		~XMLInputStream();

		bool open(const Path &path) override;
		bool object_begin(const std::string &name) override;
		bool object_end() override;

		void read(double &val) override;
		void read(int &val) override;
		void read(SizeType &val) override;
		void read(std::string &val) override;

		void read(const std::string &key, double &val) override;
		void read(const std::string &key, int &val) override;
		void read(const std::string &key, SizeType &val) override;
		void read(const std::string &key, std::string &val) override;


		bool good() const override;


		void start() override;
		void start(const std::string &name) override;

		std::string name() override;
		bool good() override;
		bool next() override;
		void finish() override;

	private:

		class Impl;
		std::unique_ptr<Impl> impl_;
	};
}

#endif //UTOPIA_XML_STREAM_HPP

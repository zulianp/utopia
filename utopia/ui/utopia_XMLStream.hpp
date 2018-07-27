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

		void read(double &val) override;
		void read(int &val) override;
		void read(SizeType &val) override;
		void read(std::string &val) override;
		void read(Serializable &val) override;
		void read(std::function<void(InputStream &)> lambda) override;

		void read(const std::string &key, double &val) override;
		void read(const std::string &key, int &val) override;
		void read(const std::string &key, SizeType &val) override;
		void read(const std::string &key, std::string &val) override;
		void read(const std::string &key, Serializable &val) override;
		void read(const std::string &key, std::function<void(InputStream &)> lambda) override;

		bool good() const override;

		SizeType size() const override;

	private:
		class Impl;
		std::unique_ptr<Impl> impl_;

		bool object_begin(const std::string &name);
		bool object_end();

		void next()   override;
		void array_start()  override;
		void array_finish() override;
	};
}

#endif //UTOPIA_XML_STREAM_HPP

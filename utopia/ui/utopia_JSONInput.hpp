#ifndef UTOPIA_JSON_STREAM_HPP
#define UTOPIA_JSON_STREAM_HPP

#include <memory>
#include "utopia_Base.hpp"
#include "utopia_Path.hpp"
#include "utopia_Input.hpp"

namespace utopia {

	class JSONInput final : public Input {
	public:
		JSONInput();
		~JSONInput();

		bool open(const Path &path);

		SizeType size() const override;
		void get(std::vector<std::shared_ptr<IConvertible>> &values) override;
		void get_all(std::function<void(Input &)> lambda) override;
		void get(const std::string &key, bool &val) override;
		void get(const std::string &key, double &val) override;
		void get(const std::string &key, int &val) override;
		void get(const std::string &key, long &val) override;
		void get(const std::string &key, unsigned long &val) override;
		// void get(const std::string &key, SizeType &val) override;
		void get(const std::string &key, std::string &val) override;
		void get(const std::string &key, Configurable &val) override;
		void get(const std::string &key, std::function<void(Input &)> lambda) override;
		bool good() const override;

	private:

		class Impl;
		std::unique_ptr<Impl> impl_;
	};
}

#endif //UTOPIA_JSON_STREAM_HPP

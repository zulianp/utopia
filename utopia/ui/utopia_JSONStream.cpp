#include "utopia_JSONStream.hpp"

#include <stack>

#include "rapidjson.h"
#include "prettywriter.h"	// for stringify JSON
#include "filestream.h"	// wrapper of C stream for prettywriter as output
#include "rapidjson/document.h"		// rapidjson's DOM-style API

#include "utopia_make_unique.hpp"
#include "utopia_Path.hpp"

#include <fstream>

using namespace rapidjson;

namespace utopia {
	class JSONInputStream::Impl {
	public:
		Impl(const Path &path)
		: current_node(nullptr), n_invalid_subtrees_(0)
		{
			std::ifstream file(path.c_str());

			if(file.good()) {

				std::stringstream ss;

				while (file.good()) {
					std::string str;
					getline(file, str);
					ss << str << "\n";
				}

				file.close();
				std::string jsonString = ss.str();

				if(doc.ParseInsitu<0>(&jsonString[0]).HasParseError())  {
					std::cerr << "Settings::fromJSON(...) : JSON Syntax error" << std::endl;
					std::cerr << doc.GetParseError() << std::endl;
				} else {
					current_node = &doc;
				}
			}
		}

		bool valid() const {

			return current_node;
		}

		bool is_invalid_subtree()
		{
			return n_invalid_subtrees_ > 0;
		}

		rapidjson::Document doc;
		rapidjson::Value *current_node;
		SizeType n_invalid_subtrees_;
	};

	bool JSONInputStream::open(const Path &path)
	{
		impl_ = make_unique<Impl>(path);
		if(!impl_->valid()) {
			impl_ = nullptr;
			return false;
		}

		return true;
	}

	JSONInputStream::~JSONInputStream() {}

	JSONInputStream::JSONInputStream() {}

	bool JSONInputStream::object_begin(const std::string &name)
	{
		//TODO
		return false;
	}

	bool JSONInputStream::object_end()
	{
		//TODO
		return false;
	}

	void JSONInputStream::read(double &val)
	{
		if(impl_->is_invalid_subtree()) return;

		//TODO
	}

	void JSONInputStream::read(int &val)
	{
		if(impl_->is_invalid_subtree()) return;
		//TODO
	}

	void JSONInputStream::read(SizeType &val)
	{
		if(impl_->is_invalid_subtree()) return;

		//TODO
	}

	void JSONInputStream::read(std::string &val)
	{
		if(impl_->is_invalid_subtree()) return;

		//TODO
	}

	void JSONInputStream::read(const std::string &key, double &val)
	{
		object_begin(key);
		read(val);
		object_end();
	}

	void JSONInputStream::read(const std::string &key, int &val)
	{
		object_begin(key);
		read(val);
		object_end();
	}

	void JSONInputStream::read(const std::string &key, SizeType &val)
	{
		object_begin(key);
		read(val);
		object_end();
	}

	void JSONInputStream::read(const std::string &key, std::string &val)
	{
		object_begin(key);
		read(val);
		object_end();
	}

	bool JSONInputStream::good() const
	{
		return impl_.get();
	}

	void JSONInputStream::start()
	{
		//TODO
	}

 	void JSONInputStream::start(const std::string &name)
 	{
 	//TODO
 	}

	std::string JSONInputStream::name()
	{
		if(impl_->is_invalid_subtree()) return "";

		//TODO
		return "";
	}

	bool JSONInputStream::good()
	{
		return !(impl_->is_invalid_subtree());
	}

	bool JSONInputStream::next()
	{
		if(impl_->is_invalid_subtree()) return false;
		//TODO
		return good();
	}

	void JSONInputStream::finish()
	{
		//TODO
		if(!good()) {
			impl_->n_invalid_subtrees_--;
		}

		object_end();
	}
}

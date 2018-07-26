#include "utopia_XMLStream.hpp"

#include <stack>

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "rapidxml_print.hpp"

#include "utopia_make_unique.hpp"
#include "utopia_Path.hpp"

using namespace rapidxml;

namespace utopia {
	class XMLInputStream::Impl {
	public:
		Impl(const Path &path)
		: f(path.c_str()), current_node(nullptr), n_invalid_subtrees_(0)
		{
			doc.parse<0>(f.data());
			current_node = &doc;
		}

		bool valid() const {
			//FIXME
			return true;
		}

		bool object_begin(const std::string &name)
		{
			if(n_invalid_subtrees_ > 0) {
				++n_invalid_subtrees_;
				return false;
			}

			// if(!current_node) {
			// 	current_node = doc.first_node(name.c_str());
			// } else {
			auto temp = current_node->first_node(name.c_str());

			if(temp) {
				current_node = temp;
			} else {
				++n_invalid_subtrees_;
			}
			// }

			return current_node;
		}

		bool is_invalid_subtree()
		{
			return n_invalid_subtrees_ > 0;
		}

		bool object_end()
		{
			if(!current_node) {
				assert(false);
				return false;
			}

			if(n_invalid_subtrees_ == 0) {
				return current_node = current_node->parent();
			} else {
				--n_invalid_subtrees_;
			}

			return false;
		}

		file<> f;
		xml_document<> doc;
		xml_node<> *current_node;
		SizeType n_invalid_subtrees_;
	};

	bool XMLInputStream::open(const Path &path)
	{
		impl_ = make_unique<Impl>(path);
		return impl_->valid();
	}

	XMLInputStream::~XMLInputStream() {}

	XMLInputStream::XMLInputStream() {}

	bool XMLInputStream::object_begin(const std::string &name)
	{
		return impl_->object_begin(name);
	}

	bool XMLInputStream::object_end()
	{
		return impl_->object_end();
	}

	void XMLInputStream::read(double &val)
	{
		if(impl_->is_invalid_subtree()) return;

		if(impl_->current_node) {
			val = atof(impl_->current_node->value());
		}
	}

	void XMLInputStream::read(int &val)
	{
		if(impl_->is_invalid_subtree()) return;

		if(impl_->current_node) {
			val = atoi(impl_->current_node->value());
		}
	}

	void XMLInputStream::read(SizeType &val)
	{
		if(impl_->is_invalid_subtree()) return;

		if(impl_->current_node) {
			val = atoi(impl_->current_node->value());
		}
	}

	void XMLInputStream::read(std::string &val)
	{
		if(impl_->is_invalid_subtree()) return;

		if(impl_->current_node) {
			val = impl_->current_node->value();
		}
	}

	void XMLInputStream::read(const std::string &key, double &val)
	{
		impl_->object_begin(key);
		read(val);
		impl_->object_end();
	}

	void XMLInputStream::read(const std::string &key, int &val)
	{
		impl_->object_begin(key);
		read(val);
		impl_->object_end();
	}

	void XMLInputStream::read(const std::string &key, SizeType &val)
	{
		impl_->object_begin(key);
		read(val);
		impl_->object_end();
	}

	void XMLInputStream::read(const std::string &key, std::string &val)
	{
		impl_->object_begin(key);
		read(val);
		impl_->object_end();
	}

	bool XMLInputStream::good() const
	{
		return impl_.get();
	}



	void XMLInputStream::start()
	{
		if(impl_->is_invalid_subtree()) {
			impl_->n_invalid_subtrees_++;
			return;
		}

		auto temp = impl_->current_node->first_node();
		if(temp) {
			impl_->current_node = temp;
		} else {
			impl_->n_invalid_subtrees_++;
		}
	}

 	void XMLInputStream::start(const std::string &name)
 	{
 		if(impl_->is_invalid_subtree()) {
 			impl_->n_invalid_subtrees_++;
 			return;
 		}

 		auto temp = impl_->current_node->first_node(name.c_str());
 		if(temp) {
 			impl_->current_node = temp;
 		} else {
 			impl_->n_invalid_subtrees_++;
 		}
 	}

	std::string XMLInputStream::name()
	{
		if(impl_->is_invalid_subtree()) return "";

		return impl_->current_node->name();
	}

	bool XMLInputStream::good()
	{
		return !(impl_->is_invalid_subtree());
	}

	bool XMLInputStream::next()
	{
		if(impl_->is_invalid_subtree()) return false;

		auto temp = impl_->current_node->next_sibling();
		if(!temp) {
			impl_->n_invalid_subtrees_++;
		} else {
			impl_->current_node = temp;
		}

		return good();
	}

	void XMLInputStream::finish()
	{
		if(!good()) {
			impl_->n_invalid_subtrees_--;
		}

		impl_->object_end();
	}
}

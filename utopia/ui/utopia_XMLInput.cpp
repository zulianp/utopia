#include "utopia_XMLInput.hpp"

#include <stack>

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "rapidxml_print.hpp"

#include "utopia_make_unique.hpp"
#include "utopia_Path.hpp"

using namespace rapidxml;

namespace utopia {
    class XMLInput::Impl {
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

            // std::cout << "object_begin: " << name <<  std::endl;

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

            // std::cout << "object_end: " << current_node->name() <<  std::endl;

            if(n_invalid_subtrees_ == 0) {
                current_node = current_node->parent();
                return true;
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

    bool XMLInput::open(const Path &path)
    {
        impl_ = make_unique<Impl>(path);
        return impl_->valid();
    }

    XMLInput::~XMLInput() {}

    XMLInput::XMLInput() {}

    bool XMLInput::object_begin(const std::string &name)
    {
        return impl_->object_begin(name);
    }

    bool XMLInput::object_end()
    {
        return impl_->object_end();
    }

    void XMLInput::get(double &val)
    {
        if(impl_->is_invalid_subtree()) return;

        if(impl_->current_node) {
            val = atof(impl_->current_node->value());
        }
    }

    void XMLInput::get(int &val)
    {
        if(impl_->is_invalid_subtree()) return;

        if(impl_->current_node) {
            val = atoi(impl_->current_node->value());
        }
    }

    void XMLInput::get(long long int &val)
    {
    	if(impl_->is_invalid_subtree()) return;

    	if(impl_->current_node) {
    		val = atoll(impl_->current_node->value());
    	}
    }

    void XMLInput::get(std::string &val)
    {
        if(impl_->is_invalid_subtree()) return;

        if(impl_->current_node) {
            val = impl_->current_node->value();
        }
    }

    void XMLInput::get(bool &val)
    {
        if(impl_->is_invalid_subtree()) return;
        static const std::string true_val = "true";

        if(impl_->current_node) {
            val = (impl_->current_node->value() == true_val);
        }
    }

    void XMLInput::get(Configurable &val)
    {
        if(impl_->is_invalid_subtree()) return;

        if(impl_->current_node) {
            val.read(*this);
        }
    }

    void XMLInput::get(const std::string &key, double &val)
    {
        impl_->object_begin(key);
        get(val);
        impl_->object_end();
    }

    void XMLInput::get(const std::string &key, int &val)
    {
        impl_->object_begin(key);
        get(val);
        impl_->object_end();
    }

    // void XMLInput::get(const std::string &key, SizeType &val)
    // {
    // 	impl_->object_begin(key);
    // 	get(val);
    // 	impl_->object_end();
    // }

    void XMLInput::get(const std::string &key, std::string &val)
    {
        impl_->object_begin(key);
        get(val);
        impl_->object_end();
    }

    void XMLInput::get(const std::string &key, bool &val)
    {
        impl_->object_begin(key);
        get(val);
        impl_->object_end();
    }

    void XMLInput::get(long &val)
    {
        if(impl_->is_invalid_subtree()) return;

        if(impl_->current_node) {
            val = atol(impl_->current_node->value());
        }
    }

    // void XMLInput::get(long long &val)
    // {
    //     if(impl_->is_invalid_subtree()) return;

    //     if(impl_->current_node) {
    //         val = atol(impl_->current_node->value());
    //     }
    // }

    void XMLInput::get(unsigned long &val)
    {
        if(impl_->is_invalid_subtree()) return;

        if(impl_->current_node) {
            val = atol(impl_->current_node->value());
        }
    }


    void XMLInput::get(const std::string &key, long &val)
    {
        impl_->object_begin(key);
        get(val);
        impl_->object_end();
    }

    // void XMLInput::get(const std::string &key, long long &val)
    // {
    //     impl_->object_begin(key);
    //     get(val);
    //     impl_->object_end();
    // }

    void XMLInput::get(const std::string &key, unsigned long &val)
    {
        impl_->object_begin(key);
        get(val);
        impl_->object_end();
    }

    void XMLInput::get(const std::string &key, long long int &val)
    {
        impl_->object_begin(key);
        get(val);
        impl_->object_end();
    }


    void XMLInput::get(const std::string &key, Configurable &val)
    {
        impl_->object_begin(key);
        get(val);
        impl_->object_end();
    }

    void XMLInput::get(std::function<void(Input &)> lambda)
    {
        if(impl_->is_invalid_subtree()) return;

        lambda(*this);
    }

    void XMLInput::get(const std::string &key, std::function<void(Input &)> lambda)
    {
        if(impl_->is_invalid_subtree()) return;

        impl_->object_begin(key);

        if(!impl_->is_invalid_subtree()) {
            lambda(*this);
        }

        impl_->object_end();
    }

    bool XMLInput::good() const
    {
        return impl_.get();
    }

    SizeType XMLInput::size() const
    {
        SizeType ret = 0;
        for(auto n = impl_->current_node->first_node(); n; n = n->next_sibling()) {
            ++ret;
        }

        return ret;
    }

    void XMLInput::array_start()
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

    void XMLInput::next()
    {
        if(impl_->is_invalid_subtree()) return;

        auto temp = impl_->current_node->next_sibling();
        if(!temp) {
            impl_->n_invalid_subtrees_++;
        } else {
            impl_->current_node = temp;
        }
    }

    void XMLInput::array_finish()
    {
        if(impl_->is_invalid_subtree()) {
            impl_->n_invalid_subtrees_--;
        }

        impl_->object_end();
    }

    void XMLInput::get_all(std::function<void(Input &)> lambda) {
        auto n = size();

        if(n == 0) return;

        array_start();

        for(SizeType i = 0; i < n; ++i) {
            lambda(*this);
            next();
        }

        array_finish();
    }
}

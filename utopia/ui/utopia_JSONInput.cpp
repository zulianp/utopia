// #include "utopia_JSONInput.hpp"

// #include <stack>

// #include "json.hpp"

// #include "utopia_make_unique.hpp"
// #include "utopia_Path.hpp"

// #include <fstream>

// namespace utopia {
// 	class JSONInput::Impl {
// 	public:
// 		using json = nlohmann::json;

// 		Impl(const Path &path)
// 		: n_invalid_subtrees_(0)
// 		{
// 			std::ifstream file(path.c_str());

// 			if(file.good()) {
// 				file >> j_;
// 				current_json_ = j_;
// 			}

// 			file.close();
// 		}

// 		bool valid() const {
// 			return !j_.empty();
// 		}

// 		bool is_invalid_subtree()
// 		{
// 			return n_invalid_subtrees_ > 0;
// 		}

// 		json j_;
// 		json current_json_;
// 		SizeType n_invalid_subtrees_;
// 	};

// 	bool JSONInput::open(const Path &path)
// 	{
// 		impl_ = make_unique<Impl>(path);
// 		if(!impl_->valid()) {
// 			impl_ = nullptr;
// 			return false;
// 		}

// 		return true;
// 	}

// 	JSONInput::~JSONInput() {}

// 	JSONInput::JSONInput() {}

// 	bool JSONInput::object_begin(const std::string &name)
// 	{
// 		//TODO


// 		return false;
// 	}

// 	bool JSONInput::object_end()
// 	{
// 		//TODO
// 		return false;
// 	}

// 	void JSONInput::read(double &val)
// 	{
// 		if(impl_->is_invalid_subtree()) return;

// 		//TODO
// 	}

// 	void JSONInput::read(int &val)
// 	{
// 		if(impl_->is_invalid_subtree()) return;
// 		//TODO
// 	}

// 	void JSONInput::read(SizeType &val)
// 	{
// 		if(impl_->is_invalid_subtree()) return;

// 		//TODO
// 	}

// 	void JSONInput::read(std::string &val)
// 	{
// 		if(impl_->is_invalid_subtree()) return;

// 		//TODO
// 	}

// 	void JSONInput::read(const std::string &key, double &val)
// 	{
// 		object_begin(key);
// 		read(val);
// 		object_end();
// 	}

// 	void JSONInput::read(const std::string &key, int &val)
// 	{
// 		object_begin(key);
// 		read(val);
// 		object_end();
// 	}

// 	void JSONInput::read(const std::string &key, SizeType &val)
// 	{
// 		object_begin(key);
// 		read(val);
// 		object_end();
// 	}

// 	void JSONInput::read(const std::string &key, std::string &val)
// 	{
// 		object_begin(key);
// 		read(val);
// 		object_end();
// 	}

// 	bool JSONInput::good() const
// 	{
// 		return impl_.get();
// 	}

// 	void JSONInput::start()
// 	{
// 		//TODO
// 	}

//  	void JSONInput::start(const std::string &name)
//  	{
//  	//TODO
//  	}

// 	std::string JSONInput::name()
// 	{
// 		if(impl_->is_invalid_subtree()) return "";

// 		//TODO
// 		return "";
// 	}

// 	bool JSONInput::good()
// 	{
// 		return !(impl_->is_invalid_subtree());
// 	}

// 	bool JSONInput::next()
// 	{
// 		if(impl_->is_invalid_subtree()) return false;
// 		//TODO
// 		return good();
// 	}

// 	void JSONInput::finish()
// 	{
// 		//TODO
// 		if(!good()) {
// 			impl_->n_invalid_subtrees_--;
// 		}

// 		object_end();
// 	}
// }

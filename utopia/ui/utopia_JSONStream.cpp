// #include "utopia_JSONStream.hpp"

// #include <stack>

// #include "json.hpp"

// #include "utopia_make_unique.hpp"
// #include "utopia_Path.hpp"

// #include <fstream>

// namespace utopia {
// 	class JSONInputStream::Impl {
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

// 	bool JSONInputStream::open(const Path &path)
// 	{
// 		impl_ = make_unique<Impl>(path);
// 		if(!impl_->valid()) {
// 			impl_ = nullptr;
// 			return false;
// 		}

// 		return true;
// 	}

// 	JSONInputStream::~JSONInputStream() {}

// 	JSONInputStream::JSONInputStream() {}

// 	bool JSONInputStream::object_begin(const std::string &name)
// 	{
// 		//TODO


// 		return false;
// 	}

// 	bool JSONInputStream::object_end()
// 	{
// 		//TODO
// 		return false;
// 	}

// 	void JSONInputStream::read(double &val)
// 	{
// 		if(impl_->is_invalid_subtree()) return;

// 		//TODO
// 	}

// 	void JSONInputStream::read(int &val)
// 	{
// 		if(impl_->is_invalid_subtree()) return;
// 		//TODO
// 	}

// 	void JSONInputStream::read(SizeType &val)
// 	{
// 		if(impl_->is_invalid_subtree()) return;

// 		//TODO
// 	}

// 	void JSONInputStream::read(std::string &val)
// 	{
// 		if(impl_->is_invalid_subtree()) return;

// 		//TODO
// 	}

// 	void JSONInputStream::read(const std::string &key, double &val)
// 	{
// 		object_begin(key);
// 		read(val);
// 		object_end();
// 	}

// 	void JSONInputStream::read(const std::string &key, int &val)
// 	{
// 		object_begin(key);
// 		read(val);
// 		object_end();
// 	}

// 	void JSONInputStream::read(const std::string &key, SizeType &val)
// 	{
// 		object_begin(key);
// 		read(val);
// 		object_end();
// 	}

// 	void JSONInputStream::read(const std::string &key, std::string &val)
// 	{
// 		object_begin(key);
// 		read(val);
// 		object_end();
// 	}

// 	bool JSONInputStream::good() const
// 	{
// 		return impl_.get();
// 	}

// 	void JSONInputStream::start()
// 	{
// 		//TODO
// 	}

//  	void JSONInputStream::start(const std::string &name)
//  	{
//  	//TODO
//  	}

// 	std::string JSONInputStream::name()
// 	{
// 		if(impl_->is_invalid_subtree()) return "";

// 		//TODO
// 		return "";
// 	}

// 	bool JSONInputStream::good()
// 	{
// 		return !(impl_->is_invalid_subtree());
// 	}

// 	bool JSONInputStream::next()
// 	{
// 		if(impl_->is_invalid_subtree()) return false;
// 		//TODO
// 		return good();
// 	}

// 	void JSONInputStream::finish()
// 	{
// 		//TODO
// 		if(!good()) {
// 			impl_->n_invalid_subtrees_--;
// 		}

// 		object_end();
// 	}
// }

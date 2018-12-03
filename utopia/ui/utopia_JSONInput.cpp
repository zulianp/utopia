#include "utopia_JSONInput.hpp"

#ifdef WITH_JSON
#include "json.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_Path.hpp"
#include "utopia_Utils.hpp"
#include <fstream>

namespace utopia {
	class JSONInput::Impl {
	public:
		Impl()
		{}

		void init(const nlohmann::json &j)
		{
			init(utopia::make_ref(j));
		}

		void init(const std::shared_ptr<const nlohmann::json> &j)
		{
			this->j = j;
		}

		Impl(const std::shared_ptr<const nlohmann::json> &j)
		: j(j)
		{}

		Impl(const Path &path)
		{
			std::ifstream file(path.c_str());
			auto j_non_const = std::make_shared<nlohmann::json>();

			if(file.good()) {
				file >> (*j_non_const);
			}

			j = std::move(j_non_const);

			file.close();
		}

		const nlohmann::json &json() const
		{
			assert(j);
			return *j;
		}

		bool good() const {
			return j && !j->empty();
		}

	private:
		std::shared_ptr<const nlohmann::json> j;

		// struct UtopiaJSONType {
		// 	/// type for (signed) integers
		// 	using number_integer_t = long;
		// 	/// type for unsigned integers
		// 	using number_unsigned_t = unsigned long;
		// 	/// type for floating-point numbers
		// 	using number_float_t = double;
		// 	/// type for strings
		// 	using string_t = std::string;
		// };

		// template<typename T>
		// class SAXJSONParser : nlohmann::json_sax<UtopiaJSONType> {
		// public:

		// 	SAXJSONParser(const T &val) 
		// 	: val_(val)
		// 	{}

		// 	// called when null is parsed
		// 	bool null()
		// 	{
		// 		return true;
		// 	}


		// 	// called when a boolean is parsed; value is passed
		// 	bool boolean(bool val)
		// 	{
		// 		val_.set(val);
		// 		return true;
		// 	}


		// 	// called when a signed or unsigned integer number is parsed; value is passed
		// 	bool number_integer(long val)
		// 	{
		// 		val_.set(val);
		// 		return true;
		// 	}

		// 	bool number_unsigned(unsigned long val)
		// 	{
		// 		val_.set(static_cast<long>(val)); //FIXME
		// 		return true;
		// 	}


		// 	// called when a floating-point number is parsed; value and original string is passed
		// 	bool number_float(double val, const std::string& s)
		// 	{
		// 		val_.set(val);
		// 		return true;
		// 	}


		// 	// called when a string is parsed; value is passed and can be safely moved away
		// 	bool string(std::string& val)
		// 	{
		// 		val_.set(val);
		// 		return true;
		// 	}


		// 	// called when an object or array begins or ends, resp. The number of elements is passed (or -1 if not known)
		// 	bool start_object(std::size_t elements)
		// 	{
		// 		return true;
		// 	}

		// 	bool end_object()
		// 	{
		// 		return true;
		// 	}

		// 	bool start_array(std::size_t elements)
		// 	{
		// 		return true;
		// 	}

		// 	bool end_array()
		// 	{
		// 		return true;
		// 	}

		// 	// called when an object key is parsed; value is passed and can be safely moved away
		// 	bool key(std::string& val)
		// 	{
		// 		// val_.set(val);
		// 		return true;
		// 	}


		// 	// called when a parse error occurs; byte position, the last token, and an exception is passed
		// 	bool parse_error(std::size_t position, const std::string& last_token, const nlohmann::detail::exception& ex)
		// 	{
		// 		assert(false);
		// 		return true;
		// 	}

		// 	inline T get() const
		// 	{	
		// 		T ret;
		// 		val_.get(ret);
		// 		return ret;
		// 	}

		// 	Convertible<T> val_;
		// };

	public:

		template<typename T>
		void get(const std::string &key, T &value) const
		{
			auto it = json().find(key);
			if(it == json().end()) {
				return;
			}

			// SAXJSONParser<T> parser(value);
			// nlohmann::json::sax_parse(*it, &parser);
			// value = parser.get();

			const auto &j = *it;
			
			if(!j.is_null()) {
				Convertible<T> c(value);
				if(j.is_boolean()) {
					c.set(j.get<bool>());

				} else if(j.is_number()) {
					
					if(c.is_double()) {
						c.set(j.get<double>());
					} else if(c.is_float()) { 
						c.set(j.get<float>());
					} else if(c.is_int()) {
						c.set(j.get<int>());
					} else if(c.is_long()) {
						c.set(j.get<long>());
					} else if(c.is_ulong()) {
						c.set(j.get<unsigned long>());
					} else if(c.is_string()) {
						c.set(j.get<double>());
					}

				} else if(j.is_string()) {
					c.set(j.get<std::string>());
				} 

				// else if(j.is_array()) {
				// 	c.set(j.get<bool>());
				// } 

				// else if(j.is_object()) {
				// 	c.set(j.get<bool>());
				// } 

				c.get(value);
			}
		}
	};

	JSONInput::JSONInput()
	{}

	JSONInput::~JSONInput()
	{}

	bool JSONInput::open(const Path &path)
	{
		impl_ = utopia::make_unique<Impl>(path);
		return impl_->good();
	}

	SizeType JSONInput::size() const
	{
		return impl_->json().size();
	}

	void JSONInput::get(std::vector<std::shared_ptr<IConvertible>> &values)
	{
		for(auto it = impl_->json().begin(); it != impl_->json().end(); ++it) {
			if(it->is_object()) continue; //FIXME?
			values.push_back(std::make_shared<Convertible<std::string>>(it->get<std::string>()));
		}
	}

	void JSONInput::get_all(std::function<void(Input &)> lambda)
	{
		for(auto it = impl_->json().begin(); it != impl_->json().end(); ++it) {
			if(!it->is_object()) continue; //FIXME?

			JSONInput child;
			child.impl_ = utopia::make_unique<Impl>(make_ref(*it));
			lambda(child);
		}
	}

	void JSONInput::get(const std::string &key, bool &val)
	{
		assert(impl_);
		impl_->get(key, val);
	}

	void JSONInput::get(const std::string &key, double &val)
	{
		assert(impl_);
		impl_->get(key, val);
	}

	void JSONInput::get(const std::string &key, int &val)
	{
		assert(impl_);
		impl_->get(key, val);
	}

	// void JSONInput::get(const std::string &key, SizeType &val)
	// {
	// 	assert(impl_);
	// 	impl_->get(key, val);
	// }

	void JSONInput::get(const std::string &key, std::string &val)
	{
		assert(impl_);
		impl_->get(key, val);
	}

	void JSONInput::get(const std::string &key, long &val)
	{
		assert(impl_);
		impl_->get(key, val);
	}

	void JSONInput::get(const std::string &key, unsigned long &val)
	{
		assert(impl_);
		impl_->get(key, val);
	}

	void JSONInput::get(const std::string &key, Configurable &val)
	{
		assert(impl_);
		auto it = impl_->json().find(key);
		if(it == impl_->json().end()) {
			// std::cout << "[Warning] key " << key << " not found in " << impl_->json().dump() << std::endl;
			return;
		}

		const auto &j = *it;

		if(!j.is_null()) {
			JSONInput child;
			child.impl_ = utopia::make_unique<Impl>(make_ref(j));
			val.read(child);
		} 
		
		// else {
		// 	std::cout << "[Warning] key " << key << " not found in " << impl_->json().dump() << std::endl;
		// }
	}

	void JSONInput::get(const std::string &key, std::function<void(Input &)> lambda)
	{
		assert(impl_);
		auto it = impl_->json().find(key);
		if(it == impl_->json().end()) {
			// std::cout << "[Warning] key " << key << " not found in " << impl_->json().dump() << std::endl;
			return;
		}

		const auto &j = *it;

		if(!j.is_null()) {
			JSONInput child;
			child.impl_ = utopia::make_unique<Impl>(make_ref(j));
			lambda(child);
		} 

		// else {
		// 	std::cout << "[Warning] key " << key << " not found in " << impl_->json().dump() << std::endl;
		// }
	}

	bool JSONInput::good() const
	{
		return (impl_) && impl_->good();
	}

}

#endif //WITH_JSON


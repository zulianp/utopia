#include "utopia_JSONInput.hpp"

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
		val = impl_->json()[key];
	}

	void JSONInput::get(const std::string &key, double &val)
	{
		assert(impl_);
		val = impl_->json()[key];
	}

	void JSONInput::get(const std::string &key, int &val)
	{
		assert(impl_);
		val = impl_->json()[key];

	}

	void JSONInput::get(const std::string &key, SizeType &val)
	{
		assert(impl_);
		val = impl_->json()[key];
	}

	void JSONInput::get(const std::string &key, std::string &val)
	{
		assert(impl_);
		val = impl_->json()[key];
	}

	void JSONInput::get(const std::string &key, Configurable &val)
	{
		assert(impl_);
		JSONInput child;
		child.impl_ = utopia::make_unique<Impl>(make_ref(impl_->json()[key]));
		val.read(child);
	}

	void JSONInput::get(const std::string &key, std::function<void(Input &)> lambda)
	{
		assert(impl_);
		JSONInput child;
		child.impl_ = utopia::make_unique<Impl>(make_ref(impl_->json()[key]));
		lambda(child);
	}

	bool JSONInput::good() const
	{
		return (impl_) && impl_->good();
	}


}

#ifndef UTOPIA_INPUT_PARAMETERS_HPP
#define UTOPIA_INPUT_PARAMETERS_HPP

#include "utopia_Input.hpp"
#include "utopia_Convertible.hpp"
#include "utopia_make_unique.hpp"

#include <map>
#include <cassert>
#include <string>
#include <memory>

namespace utopia {
	class InputParameters /*final*/ : public Input {
	public:
		inline bool empty() const
		{
			return nodes_.empty() && values_.empty();
		}
		
		inline SizeType size() const override
		{
			return nodes_.size() + values_.size();
		}

		inline void read(const std::string &key, bool &val) override
		{
			aux_read(key, val);
		}

		inline void read(const std::string &key, double &val) override
		{
			aux_read(key, val);
		}

		inline void read(const std::string &key, int &val) override
		{
			aux_read(key, val);
		}

		inline void read(const std::string &key, SizeType &val) override
		{
			aux_read(key, val);
		}

		inline void read(const std::string &key, std::string &val) override
		{
			aux_read(key, val);
		}

		inline void read(const std::string &key, Configurable &val) override
		{
			auto node_ptr = node(key);

			if(node_ptr) {
				val.read(*node_ptr);
			}
		}

		void read_all(std::function<void(Input &)> lambda) override
		{
			for(auto n : nodes_) {
				lambda(*n.second);
			}
		}

		inline bool good() const override { return true; }

		std::shared_ptr<Input> node(const std::string &key) const
		{
			auto it = nodes_.find(key);
			if(it != nodes_.end()) {
				return it->second;
			}

			assert(false);
			return nullptr;
		}

		inline void set(const std::string &key, const bool &val)
		{
			aux_set(key, val);
		}

		inline void set(const std::string &key, const double &val)
		{
			aux_set(key, val);
		}

		inline void set(const std::string &key, const int &val)
		{
			aux_set(key, val);
		}

		inline void set(const std::string &key, const SizeType &val)
		{
			aux_set(key, val);
		}

		inline void set(const std::string &key, const std::string &val)
		{
			aux_set(key, val);
		}

		inline void set(const std::string &key, std::shared_ptr<Input> &in)
		{
			nodes_[key] = in;
		}

	private:
		std::map<std::string, std::unique_ptr<IConvertible>> values_;
		std::map<std::string, std::shared_ptr<Input>> nodes_;

		template<typename Out>
		void aux_read(const std::string &key, Out &out) const {
			auto it = values_.find(key);

			if(it != values_.end()) {
				it->second->get(out);
			} //else do nothing
		}

		template<typename In>
		void aux_set(const std::string &key, const In &in) {
			values_[key] = make_unique<Convertible<In>>(in);
		}
	};
}


#endif //UTOPIA_INPUT_PARAMETERS_HPP

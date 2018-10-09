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
		inline void read(const std::string &key, bool &val) 
		{
			aux_read(key, val);
		}

		inline void read(const std::string &key, double &val) 
		{
			aux_read(key, val);
		}

		inline void read(const std::string &key, int &val) 
		{
			aux_read(key, val);
		}

		inline void read(const std::string &key, SizeType &val) 
		{
			aux_read(key, val);
		}

		inline void read(const std::string &key, std::string &val) 
		{
			aux_read(key, val);
		}

		inline void read(const std::string &key, Configurable &val) 
		{
			auto node_ptr = node(key);
			
			if(node_ptr) {
				val.read(*node_ptr);
			}
		}

		inline bool good() const { return true; }

		std::shared_ptr<Input> node(const std::string &key) const
		{
			auto it = nodes_.find(key);
			if(it != nodes_.end()) {
				return it->second;
			}

			assert(false);
			return nullptr;
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
	};
}


#endif //UTOPIA_INPUT_PARAMETERS_HPP

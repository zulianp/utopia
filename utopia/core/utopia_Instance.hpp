//
// Created by Patrick Zulian on 26/05/15.
//

#ifndef UTOPIA_UTOPIA_INSTANCE_HPP
#define UTOPIA_UTOPIA_INSTANCE_HPP

#include <map>
#include <string>

namespace utopia {
    class Utopia {
    public:
       static void Init(int argc, char *argv[]);
       static int Finalize();

       inline std::string get(const std::string &key) const
       {
       		auto it = settings_.find(key);
       		if(it == settings_.end()) return "";
       		return it->second;
       }

       inline void set(const std::string &key, const std::string &value)
       {
       		settings_[key] = value;
       }

       static Utopia &Instance();

      bool verbose() const;

   private:
   		Utopia();
       	std::map<std::string, std::string> settings_;
    };
}
#endif //UTOPIA_UTOPIA_INSTANCE_HPP

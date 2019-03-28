#ifndef UTOPIA_CODE_TEMPLATE_HPP
#define UTOPIA_CODE_TEMPLATE_HPP


#include <map>
#include <iostream>
#include <vector>

namespace utopia {
    static const char * TPL_LIST_ATTR = "[list]";

    class CodeTemplate {
    public:
        void init();
        void map(const std::string &key, const std::string &value);
        // void map_dynamic(const std::string &key, const std::string &value);
        bool parse(const std::string &tpl, std::ostream &os) const;
        bool parse(const std::string &tpl, std::string &result) const;


        CodeTemplate();
    private:

        typedef std::map<std::string, std::vector<std::string> > ListMap;

        std::string _varStart;
        std::string _varEnd;
        std::map<std::string, std::string> _mappings;
        ListMap _listmappings;

        bool is_list(const std::string &var) const;
    };
}


#endif //UTOPIA_CODE_TEMPLATE_HPP

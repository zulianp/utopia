
#include "utopia_opencl_CodeTemplate.hpp"
#include "utopia_Base.hpp"

#include <sstream>
#include <assert.h>

namespace utopia {

    void CodeTemplate::init()
    {
        _varStart = "<@";
        _varEnd = "@>";
    }

    CodeTemplate::CodeTemplate()
    {
        init();
    }

    void CodeTemplate::map(const std::string &key, const std::string &value)
    {
        if(is_list(key)) {
            _listmappings[key].push_back(value);
            return;
        }

        _mappings[key] = value;
    }


    bool CodeTemplate::is_list(const std::string &var) const
    {
        return var.find(TPL_LIST_ATTR) != std::string::npos;
    }

    bool CodeTemplate::parse(const std::string &tpl, std::string &result) const {

        std::stringstream ss;
        if(!parse(tpl, ss)) {
            return false;
        }

        result = ss.str();
        return true;
    }

    bool CodeTemplate::parse(const std::string &tpl, std::ostream &os) const {
        size_t endPos = std::string::npos;
        size_t startPos = tpl.find(_varStart);


        if(startPos == std::string::npos) {
            std::cerr << "------------------------------------------------------" << std::endl;
            std::cerr << "[Error] sintax error" << std::endl;
            std::cerr << "------------------------------------------------------" << std::endl;
            return false;
        }

        while(startPos != std::string::npos) {
            if(endPos != std::string::npos) {
                os << tpl.substr(endPos+_varEnd.size(), startPos-endPos-_varEnd.size());
            } else {
                os << tpl.substr(0, startPos);
            }

            size_t startOffset = startPos + _varStart.size();

            endPos = tpl.find(_varEnd, startOffset);
            assert(endPos != std::string::npos);

            if(endPos == std::string::npos) {
                std::cerr << "------------------------------------------------------" << std::endl;
                std::cerr << "[Error] sintax error" << std::endl;
                std::cerr << "------------------------------------------------------" << std::endl;
                return false;
            }





            const std::string var = tpl.substr(startOffset, endPos-startOffset);

            if(is_list(var)) {
                ListMap::const_iterator it = _listmappings.find(var);

                if(it == _listmappings.end()) {
                    // std::cerr << "------------------------------------------------------" << std::endl;
                    // std::cerr << "[Warning] non existing template variable (empty list): " << var << std::endl;
                    // std::cerr << "------------------------------------------------------" << std::endl;
                    // assert(false);
                    // return false;
                } else {
                    const std::vector<std::string> &list = it->second;
                    for(std::vector<std::string>::const_iterator lIt = list.begin(); lIt != list.end(); ++lIt) {
                        os << *lIt;
                        if(lIt+1 != list.end()) {
                            os << ", ";
                        }
                    }
                }

            } else {
                std::map<std::string, std::string>::const_iterator it = _mappings.find(var);

                if(it == _mappings.end()) {
                    std::cerr << "------------------------------------------------------" << std::endl;
                    std::cerr << "[Error] non existing template variable: " << var << std::endl;
                    std::cerr << "------------------------------------------------------" << std::endl;
                    return false;
                }

                os << it->second;
            }





            startPos = tpl.find(_varStart, endPos + _varEnd.size());
        }


        os << tpl.substr(endPos+_varEnd.size(), tpl.size()-endPos-_varEnd.size());
        return true;
    }
}

#include "utopia_Utils.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <ostream>
#include <iostream>

namespace utopia 
{

	bool read(const std::string &path, std::string &str)
	{
		std::stringstream content;
		std::string line;                                
		std::ifstream file (path.c_str());

		if (file.is_open())
		{
			while ( file.good() )
			{
				getline (file,line);
				content << line << "\n";
			}
			file.close();
		}
		else {

			std::cerr<<"Unable the open "<<path<<std::endl;
			return false;
		}

		str = content.str();
		return true;
	}

	bool write(const std::string &path, const std::string &str)
	{
		std::ofstream os(path.c_str());
		if(os.is_open()) {
			os << str;
			os.close();
			return true;
		} else {
			std::cerr<<"Unable the open "<<path<<std::endl;
			os.close();
			return false;
		}
	}

	void disp(const std::string &str)
	{
		if(mpi_world_rank() == 0) 	
			std::cout << str << std::endl;
	}

	void disp(const char * str)
	{
		if(mpi_world_rank() == 0) 	
			std::cout << str << std::endl;
	}


	/**
	 * @brief      Splits string according to requested delimiter
	 * @addtogroup Utilities
	 * @param[in]  s      string
	 * @param[in]  delim  The delimiter
	 *
	 * @return     Vector of splits 
	 */
    std::vector<std::string> split_string(const std::string &s, char delim) 
	{
	    std::vector<std::string> v;
	    auto i = 0;
	    auto pos = s.find(delim);
	    while (pos != std::string::npos)
	    {
	        v.push_back(s.substr(i, pos-i));
	        i = ++pos;
	        pos = s.find(delim, pos);
	        if (pos == std::string::npos)
	            v.push_back(s.substr(i, s.length()));
	    }
	    if(v.size() == 0) // if only one word is in the string
	        v.push_back(s);
	    return v;
	}


}


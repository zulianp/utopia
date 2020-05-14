/*
* @Author: kopanicakova
* @Date:   2017-07-04 00:18:50
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-02-12 16:21:50
*/
#include "utopia_Utils.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <ostream>
#include <iostream>

namespace utopia
{
    bool is_matlab_file(const std::string &path)
    {
        size_t pos = path.find_last_of(".");

        bool is_matlab = false;
        if(pos != std::string::npos) {
            is_matlab = path.substr(pos, path.size()-pos) == ".m";
        }

        return is_matlab;
    }


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


    template<>
     void CSVWriter::write_table_row(const std::vector<double> vars)
     {
        if (pFile != nullptr && mpi_world_rank() == 0) {
            for (std::vector<double>::size_type i = 0; i < vars.size(); i++) {
                if (i < vars.size() - 1)
                    fprintf(pFile, "%g,", vars[i]);
                else
                    fprintf(pFile, "%g", vars[i]);
            }
            fprintf(pFile, "\n");
        }
     }

     template<>
     void CSVWriter::write_table_row(const std::vector<int> vars)
     {
         if (pFile != nullptr && mpi_world_rank() == 0) {
             for(std::vector<double>::size_type i = 0;  i < vars.size(); i++ )
             {
                 if(i < vars.size() -1)
                     fprintf (pFile,"%d,", vars[i]);
                 else
                     fprintf (pFile,"%d", vars[i]);
             }
             fprintf (pFile, "\n");
         }
     }

     template<>
     void CSVWriter::write_table_row(const std::vector<unsigned long> vars)
     {
         if (pFile != nullptr && mpi_world_rank() == 0) {
             for(std::vector<double>::size_type i = 0;  i < vars.size(); i++ )
             {
                 if(i < vars.size() -1)
                     fprintf (pFile,"%lu,", vars[i]);
                 else
                     fprintf (pFile,"%lu", vars[i]);
             }
             fprintf (pFile, "\n");
         }
     }

     template<>
     void CSVWriter::write_table_row(const std::vector<char> vars)
     {
         if (pFile != nullptr && mpi_world_rank() == 0) {
             for(std::vector<double>::size_type i = 0;  i < vars.size(); i++ )
             {
                 if(i < vars.size() -1)
                     fprintf (pFile,"%c,", vars[i]);
                 else
                     fprintf (pFile,"%c", vars[i]);
             }
             fprintf (pFile, "\n");
         }
     }

     template<>
     void CSVWriter::write_table_row(const std::vector<std::string> vars)
     {
         if (pFile != nullptr && mpi_world_rank() == 0) {
             for(std::vector<double>::size_type i = 0;  i < vars.size(); i++ )
             {
                 if(i < vars.size() -1)
                     fprintf (pFile,"%s,", vars[i].c_str());
                 else
                     fprintf (pFile,"%s", vars[i].c_str());
             }
             fprintf (pFile, "\n");
         }
     }


}


#ifndef UTOPIA_SOLVER_PRINT_Info_HPP
#define UTOPIA_SOLVER_PRINT_Info_HPP


#include <iomanip>
#include <limits>
#include <chrono>
#include "utopia_Utils.hpp"
#include "utopia_Traits.hpp"

namespace utopia
{
    /**
     * @brief      The class helping to print-out information about solver: initialization messages, prinitning iteration status, time-stats and exit/convergance messages.
     *             It also helps pass solution and informations about solve back into FEM packages.
     */
    class PrintInfo
    {
        typedef double Scalar;

    public:

        static void print_init(const std::string & method, const std::vector<std::string> status_variables)
        {
            if(mpi_world_rank() == 0)
            {
                std::cout<<"  \n"; 
                std::cout << std::setw(10)  << std::right << std::string(80, '_') <<  std::setw(15)<< std::endl;
                std::cout<<  std::setw(45) << std::right << " utopia:: "<< method << std::setw(15)<< std::endl;
                std::cout << std::setw(10)  << std::right << std::string(80, '_') <<  std::setw(15)<< std::endl;
                std::cout << std::endl;

                for(auto item =  status_variables.begin(); item!= status_variables.end(); item++ )
                    std::cout << std::setw(27) << std::right << *item ;
                std::cout<<std::endl;

                auto n = status_variables.size();
                for(size_t i = 0; i < n; i ++ )
                    std::cout << std::setw(27) << std::right << std::string(10, '-');
                std::cout<<std::endl;
            }
        }


        static void print_iter_status(const std::vector<Scalar> status_variables)
        {
            if(mpi_world_rank() == 0)
            {
                for(auto item =  status_variables.begin(); item!= status_variables.end(); item++ )
                    std::cout << std::setw(27) <<  std::right << std::scientific << *item;
                std::cout << std::endl;
            }
        }


        template<class T>
        static void print_iter_status(const T & it, const std::vector<Scalar> scalar_vars)
        {
            if(mpi_world_rank() == 0)
            {
                std::cout.precision(15);
                std::cout << std::setw(27) << std::right << it;
                for(auto item =  scalar_vars.begin(); item!= scalar_vars.end(); item++ )
                    std::cout << std::setw(27) << std::right << std::scientific <<  *item;
                std::cout<<std::endl;
            }
        }


    private:
        PrintInfo()
        {

        }


    };

}
#endif //UTOPIA_SOLVER_PRINT_Info_HPP

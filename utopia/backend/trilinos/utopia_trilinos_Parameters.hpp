
#ifndef UTOPIA_TRILINOS_PARAMETERS_HPP
#define UTOPIA_TRILINOS_PARAMETERS_HPP

#include "Teuchos_ParameterList.hpp"
#include "utopia_Parameters.hpp"


namespace utopia
{
// creare oggetto nel back-end di rilinos che copia in un oggetto teuchos i paramteri di Parameter

/**
* @brief      This class keeps track on all parameters, that we have in linear solvers, it is derived from the utopia_Parameters class.
*
*/
class Trilinos_Parameters
    {
        typedef double Scalar;
        typedef long SizeType;

    public:
        // copy default parameters in ParameterLis obj
	Parameters() // default constructor
            {

            } //Parameters()

        Parameters(Teuchos::ParameterList& params)
            {
            belosParams = params;
            }


    protected:
	Teuchos::ParameterList belosParams;
    } //Trilinos_Parameters
}
#endif //UTOPIA_TRILINOS_PARAMETERS_HPP

#include "utopia.hpp"

#include <memory>
#include <iostream>
#include <sstream>
#include <ctime>

#include "utopia_Allocations.hpp"

int main(const int argc, char *argv[])
{
    using namespace std;
    using namespace utopia;
    Utopia::Init(argc, argv);
    //TODO add app runner here
    return Utopia::Finalize();
}

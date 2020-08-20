
#include <iostream>
#include "petscversion.h"

int main() {
  std::cout << PETSC_VERSION_MAJOR << "."
        << PETSC_VERSION_MINOR << "."
        << PETSC_VERSION_SUBMINOR;
  return 0;
}

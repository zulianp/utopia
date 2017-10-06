#include "utopia_Base.hpp"
#include "utopia_TrilinosTest.hpp"

#ifdef WITH_TRILINOS


#include "utopia.hpp"

using namespace Teuchos;

void trilinos_example_test()
{
	RCP<const Comm<int> > comm = rcp(new MpiComm<int>(MPI_COMM_WORLD));
  if (comm->getRank () == 0) {
    // On (MPI) Process 0, print out the Tpetra software version.
    std::cout << Tpetra::version () << std::endl << std::endl;
  }
}

namespace utopia {
	void run_trilinos_test()
	{
		UTOPIA_UNIT_TEST_BEGIN("TrilinosTest");
		UTOPIA_RUN_TEST(trilinos_example_test);
		UTOPIA_UNIT_TEST_END("TrilinosTest");
	}
}

#else //WITH_TRILINOS
namespace utopia {
	void run_trilinos_test() {}
}
#endif //WITH_TRILINOS
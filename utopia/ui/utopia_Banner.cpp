#include "utopia_Banner.hpp"

#include "utopia_IOStream.hpp"
#include "utopia_MPI.hpp"

static const char* utopia_logo_ =
    "--------------------------------------------------------------------------------\n"  //
    "--------------------------------------------------------------------------------\n"  //
    "                                                  00                            \n"  //
    "                                                  00                            \n"  //
    "                                                  00                            \n"  //
    "                                                  00                            \n"  //
    "000        00   0000000000000     1000001      0000000001     000       000     \n"  //
    "000        00   0000000000000   000000000000   0000000000000  000      00000    \n"  //
    "000        00        000       000        00      00    1000  000     000 000   \n"  //
    "000       100        000       000        00   00000000001    000    000   000  \n"  //
    " 000      000        000        000      000      00          000   000     000 \n"  //
    "  1000000001         000         1000000001       00          000  001       100\n"  //
    "                                                  00                            \n"  //
    "--------------------------------------------------------------------------------\n"  //
    "--------------------------------------------------------------------------------\n";

static const char* info =
    "Git-repository:\thttps://bitbucket.org/zulianp/utopia/src/master/\n"
    "License:\thttps://bitbucket.org/zulianp/utopia/src/master/LICENSE.md\n"
    "Credits\t\tIf you use utopia for your research use `--citations` option "
    "to your runs to see which papers have to be cited.\n";

namespace utopia {
    void Banner::welcome() {
        if (mpi_world_rank() == 0) {
            utopia::out() << utopia_logo_;
            utopia::out() << info;
        }
    }

    void Banner::bye() {
        if (mpi_world_rank() == 0) {
            utopia::out() << "Bye!\n";
        }
    }
}  // namespace utopia

#include "utopia_Banner.hpp"

#include "utopia_MPI.hpp"
#include "utopia_IOStream.hpp"

static const char * utopia_logo_ =
"------------------------------------------------------------------------------------------------------------------------\n" //
"------------------------------------------------------------------------------------------------------------------------\n" //
"                                                                           000                                          \n" //
"                                                                           000                                          \n" //
"                                                                           000                                          \n" //
"                                                                           000                                          \n" //
"                                                                           000                                          \n" //
"                                                                           000                                          \n" //
"0000            0000    0000000000000000000         100000001         00000000000001         0000           0000        \n" //
"0000            0000    0000000000000000000     00000000000000000     0000000000000000000    0000          000000       \n" //
"0000            0000           0000            0000           0000         000       0000    0000         00000000      \n" //
"0000            0000           0000            0000           0000         000      00000    0000        0000  0000     \n" //
"0000            0000           0000            0000           0000    000000000000000000     0000       0000    0000,   \n" //
"0000            0000           0000            0000           0000         000               0000     00000      00000  \n" //
" 00000         00000           0000            00000         00000         000               0000    00000        00000 \n" //
"  0000000000000000             0000             00000000000000000          000               0000   00001           0000\n" //
"                                                                           000                                          \n" //
"                                                                           000                                          \n" //
"                                                                           000                                          \n" //
"                                                                           000                                          \n" //
"------------------------------------------------------------------------------------------------------------------------\n" //
"------------------------------------------------------------------------------------------------------------------------\n"; 

static const char * info = 
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
}

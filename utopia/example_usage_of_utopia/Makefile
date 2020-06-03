UTOPIA_INSTALLATION_PATH=$(UTOPIA_DIR)

include $(UTOPIA_INSTALLATION_PATH)/config/utopia-config.makefile

INCLUDES:= $(UTOPIA_INCLUDES)
LIBS:=$(UTOPIA_LIBRARIES)
DEFS:=$(UTOPIA_DEFS)
CXX:=mpicc
CC:=mpic++
CXX_FLAGS:=$(UTOPIA_CXX_FLAGS) $(CXX_FLAGS) -std=c++14

all: hello

hello: ../examples/example_hello_world.cpp $(LIB_UTOPIA)
	$(CXX) $(CXX_FLAGS) $(DEFS) ../examples/example_hello_world.cpp $(INCLUDES) -o hello $(LIBS)

clean:
	rm -f hello rwt fem1D rosenbrock utopia_mult_exp

run_hello: hello
	mpiexec -np 3 ./hello

.PHONY: clean run_hello


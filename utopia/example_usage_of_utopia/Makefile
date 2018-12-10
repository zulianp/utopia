UTOPIA_INSTALLATION_PATH=$(UTOPIA_DIR)

include $(UTOPIA_INSTALLATION_PATH)/config/utopia-config.makefile

INCLUDES:= $(UTOPIA_INCLUDES)
LIBS:=$(UTOPIA_LIBRARIES)
DEFS:=$(UTOPIA_DEFS)	
CXX:=$(UTOPIA_CXX_COMPILER)	
CC:=$(UTOPIA_C_COMPILER)	
CXX_FLAGS:=$(UTOPIA_CXX_FLAGS) $(CXX_FLAGS)	-std=c++11
LIB_UTOPIA:=$(UTOPIA_INSTALLATION_PATH)/lib/libutopia.a

all: assemble rwt fem1D rosenbrock utopia_mult_exp each_loop

assemble: utopia_assembly_example.cpp $(LIB_UTOPIA)
	$(CXX) $(CXX_FLAGS) $(DEFS) utopia_assembly_example.cpp $(INCLUDES) -o assemble $(LIBS)

fem1D: utopia_solution_of_1D_FEM.cpp $(LIB_UTOPIA)
	$(CXX) $(CXX_FLAGS) $(DEFS) utopia_solution_of_1D_FEM.cpp $(INCLUDES) -o fem1D $(LIBS)

rwt: utopia_read_write_transform_example.cpp $(LIB_UTOPIA)
	$(CXX) $(CXX_FLAGS) $(DEFS) utopia_read_write_transform_example.cpp $(INCLUDES) -o rwt $(LIBS)

rosenbrock: utopia_nonlinear_solve_example.cpp $(LIB_UTOPIA)
	$(CXX) $(CXX_FLAGS) $(DEFS) utopia_nonlinear_solve_example.cpp $(INCLUDES) -o rosenbrock $(LIBS)	

utopia_mult_exp: utopia_mult_exp.cpp $(LIB_UTOPIA)
	$(CXX) $(CXX_FLAGS) -Ofast $(DEFS) utopia_mult_exp.cpp $(INCLUDES) -o utopia_mult_exp $(LIBS)	

each_loop: utopia_each_loop_example.cpp $(LIB_UTOPIA)
	$(CXX) $(CXX_FLAGS) $(DEFS) -O3 -DNDEBUG utopia_each_loop_example.cpp $(INCLUDES) -o each_loop $(LIBS)

clean:
	rm -f assemble rwt fem1D rosenbrock utopia_mult_exp

run_assemble: assemble
	mpiexec -np 3 ./assemble	

run_rwt: rwt
	mpiexec -np 3 ./rwt	

run_fem1D: fem1D
	mpiexec -np 2 ./fem1D

run_rosenbrock: rosenbrock
	mpiexec -np 1 ./rosenbrock	
	
.PHONY: clean run
	

UTOPIA_FE_INSTALLATION_PATH=../build/
include $(UTOPIA_FE_INSTALLATION_PATH)/config/utopia_fe_config.makefile

INCLUDES:=$(UTOPIA_FE_INCLUDES)
LIBS:=$(UTOPIA_FE_LIBRARIES)
DEFS:=$(UTOPIA_FE_DEFS)	
CXX:=$(UTOPIA_FE_CXX_COMPILER)	
CC:=$(UTOPIA_FE_C_COMPILER)	
CXX_FLAGS:=$(UTOPIA_FE_CXX_FLAGS) $(CXX_FLAGS)	


all: assemble 

assemble: utopia_fe_assembly_example.cpp  $(UTOPIA_FE_INSTALLATION_PATH)/lib/libutopia_fe.a
	$(CXX) $(CXX_FLAGS) $(DEFS) utopia_fe_assembly_example.cpp $(INCLUDES) -o assemble $(LIBS)

clean:
	rm -f assemble 

run_assemble: assemble
	mpiexec -np 3 ./assemble	
	
.PHONY: clean run
	
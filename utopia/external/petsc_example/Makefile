# UTOPIA_INSTALLATION_PATH=$(MOONOLITH_ROOT)/utopia/build
UTOPIA_INSTALLATION_PATH=/Users/alenakopanicakova/MOOSE/moonolith/utopia_installation

include $(UTOPIA_INSTALLATION_PATH)/config/utopia-config.makefile

INCLUDES:=-I $(UTOPIA_INSTALLATION_PATH)/include $(UTOPIA_INCLUDES)
LIBS:=-L$(UTOPIA_INSTALLATION_PATH)/lib -lutopia $(UTOPIA_LIBRARIES)
DEFS:=$(UTOPIA_DEFS)	
CXX:=$(UTOPIA_CXX_COMPILER)	
CC:=$(UTOPIA_C_COMPILER)	
CXX_FLAGS:=$(UTOPIA_CXX_FLAGS) $(CXX_FLAGS)	
LIB_UTOPIA:=$(UTOPIA_INSTALLATION_PATH)/lib/libutopia.a

all: nonlinear_ex  

nonlinear_ex: Bratu_example.cpp $(LIB_UTOPIA)
	$(CXX) $(CXX_FLAGS) $(DEFS) Bratu_example.cpp $(INCLUDES) -o nonlinear_ex $(LIBS)

clean:
	rm -f nonlinear_ex 

run_ex: nonlinear_ex
	mpiexec -np 1 ./nonlinear_ex	

run_nonlinear_ex: nonlinear_ex
	mpiexec -np 3 ./nonlinear_ex		
	
.PHONY: clean run
	

UTOPIA_FE_INSTALLATION_PATH=../build/
include $(UTOPIA_FE_INSTALLATION_PATH)/config/utopia_fe_config.makefile

INCLUDES:=$(UTOPIA_FE_INCLUDES)
LIBS:=$(UTOPIA_FE_LIBRARIES)
DEFS:=$(UTOPIA_FE_DEFS)	
CXX:=$(UTOPIA_FE_CXX_COMPILER)	
CC:=$(UTOPIA_FE_C_COMPILER)	
CXX_FLAGS:=$(UTOPIA_FE_CXX_FLAGS) $(CXX_FLAGS) -g


all: assemble  
# sphere_mg

assemble: utopia_fe_assembly_example.cpp  $(UTOPIA_FE_INSTALLATION_PATH)/lib/libutopia_fe.a
	$(CXX) $(CXX_FLAGS) $(DEFS) utopia_fe_assembly_example.cpp $(INCLUDES) -o assemble $(LIBS)

# sphere_mg: sphere_mg_test.cpp  $(UTOPIA_FE_INSTALLATION_PATH)/lib/libutopia_fe.a
# 	$(CXX) $(CXX_FLAGS) $(DEFS) sphere_mg_test.cpp $(INCLUDES) -o sphere_mg $(LIBS)

clean:
	rm -f assemble sphere_mg

run_assemble: assemble
	mpiexec -np 3 ./assemble	

# run_sphere_mg: sphere_mg
# 	mpiexec -np 3 ./sphere_mg
	
.PHONY: clean run
	

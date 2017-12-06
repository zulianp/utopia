

list(APPEND LOCAL_HEADERS 
	utopia_Biomechanics.hpp
	utopia_ContactSimParams.hpp
	utopia_ContactTest.hpp
	utopia_ContactProblem.hpp
	utopia_FEDSLBaseExamples.hpp
	utopia_FEDSLBaseSolverExamples.hpp
	utopia_FEDSLLeastSquaresExamples.hpp
	utopia_FEDSLMortarExamples.hpp
	utopia_GeometryTest.hpp
	utopia_MixedFESpaceExample.hpp
	utopia_SemigeometricMultigridTest.hpp
	utopia_TimeDiffExamples.hpp
	utopia_UtopiaFETests.hpp
	utopia_VolumeTransferBenchmark.hpp
	utopia_FormEvalTest.hpp
	utopia_FSITest.hpp
	utopia_SDCTest.hpp
	utopia_WearEstimator.hpp
	)

list(APPEND LOCAL_SOURCES 
	utopia_Biomechanics.cpp
	utopia_ContactTest.cpp
	utopia_ContactProblem.cpp
	utopia_FEDSLBaseExamples.cpp
	utopia_FEDSLBaseSolverExamples.cpp
	utopia_FEDSLLeastSquaresExamples.cpp
	utopia_FEDSLMortarExamples.cpp
	utopia_GeometryTest.cpp
	utopia_MixedFESpaceExample.cpp
	utopia_SemigeometricMultigridTest.cpp
	utopia_TimeDiffExamples.cpp
	utopia_UtopiaFETests.cpp
	utopia_VolumeTransferBenchmark.cpp
	utopia_FormEvalTest.cpp
	utopia_LibMeshbackendTest.cpp
	utopia_FSITest.cpp
	utopia_SDCTest.cpp
	utopia_WearEstimator.cpp
	)

add_library(utopia_fe_tests STATIC ${LOCAL_SOURCES})
target_link_libraries(utopia_fe_tests transfer utopia_fe_visualization utopia_fe_homemade utopia_fe_libmesh_backend)
target_include_directories(utopia_fe_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) 


install(TARGETS utopia_fe_tests
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_tests" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${UTOPIA_FE_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -g")

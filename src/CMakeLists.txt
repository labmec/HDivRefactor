# let us create a config file with the directory containing the input files
set(INPUTSDIR ${PROJECT_SOURCE_DIR}/input)
set(CONFIG_H "${CMAKE_CURRENT_SOURCE_DIR}/dirs_config.h")
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/dirs_config.h.in" ${CONFIG_H})

target_include_directories(libproject PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SUBFLOW_SRC CONFIGURE_DEPENDS "*.cpp" "*.c" "*.hpp")
file(GLOB SUBFLOW_HEADERS CONFIGURE_DEPENDS "*.h")

target_sources(libproject PRIVATE ${SUBFLOW_SRC} ${SUBFLOW_HEADERS})
cmake_minimum_required(VERSION 2.8)

project(bam-somaticsniper)

set(SOURCES 
    main.c
)

set(EXECUTABLE_NAME bam-somaticsniper)
add_executable(${EXECUTABLE_NAME} ${SOURCES})
target_link_libraries(${EXECUTABLE_NAME} sniper)
set_target_properties(${EXECUTABLE_NAME} PROPERTIES PACKAGE_OUTPUT_NAME ${EXECUTABLE_NAME}${EXE_VERSION_SUFFIX})
install(TARGETS ${EXECUTABLE_NAME} DESTINATION bin/)

# Add the source files for the tests executable
add_executable(tests
    main.c
)

target_include_directories(tests PUBLIC ../src)

# Link the shared library to the tests executable
target_link_libraries(tests gflib)

# Add compiler flags
target_compile_options(tests PRIVATE -Wall -g)

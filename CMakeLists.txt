cmake_minimum_required(VERSION 3.10)

project(ff-libc)

# Add the source files for the static library
add_library(gflib SHARED
    src/gf.c
    src/poly.c
)

# Add the source files for the tests executable
add_executable(tests
    tests/main.c
)

# Link the shared library to the tests executable
target_link_libraries(tests gflib)

# Add compiler flags
target_compile_options(tests PRIVATE -Wall -g)

enable_testing()
add_test(NAME gflib_tests COMMAND tests)

# Add Clang format
file(GLOB_RECURSE ALL_SOURCE_FILES *.c *.h)
add_custom_target(format COMMAND clang-format -style=Google -i ${ALL_SOURCE_FILES})

# Add Valgrind
add_custom_target(memcheck COMMAND valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./tests DEPENDS tests)

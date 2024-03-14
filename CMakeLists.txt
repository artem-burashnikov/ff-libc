cmake_minimum_required(VERSION 3.10)
project(ff-libc)

# Include subdirectories
add_subdirectory(src)
add_subdirectory(tests)

enable_testing()
add_test(NAME gflib_tests COMMAND tests)

# Add Clang format
file(GLOB_RECURSE ALL_SOURCE_FILES *.c *.h)
add_custom_target(format COMMAND clang-format -style=Google -i ${ALL_SOURCE_FILES})

# Add Valgrind
add_custom_target(memcheck COMMAND valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./tests/tests DEPENDS tests)


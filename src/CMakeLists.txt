# Add static library sources.
add_library(gf STATIC
    GF.c
    poly.c
    utils.c
)

target_include_directories(gf PUBLIC ${CMAKE_SOURCE_DIR}/src)

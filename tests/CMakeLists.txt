include(CTest)

# generate an image using the latest build of cleanpt
add_test(
    NAME generate_test_image
    COMMAND cleanpt 4
)

# then compare the output with a ground truth image
add_test(
    NAME generated_image
    COMMAND ${CMAKE_COMMAND} -E compare_files
        ${CMAKE_BINARY_DIR}/tests/image.ppm
        ${CMAKE_SOURCE_DIR}/tests/image4.ppm
)

# TODO add proper unit tests once the code is split between a library and an app
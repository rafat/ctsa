add_executable(acftest acftest.c)

target_link_libraries(acftest ctsalib m)

add_executable(arimatest arimatest.c)

target_link_libraries(arimatest ctsalib m)

add_executable(artest artest.c)

target_link_libraries(artest ctsalib m)

add_executable(artest2 artest2.c)

target_link_libraries(artest2 ctsalib m)

add_executable(pacftest pacftest.c)

target_link_libraries(pacftest ctsalib m)

add_executable(sarimatest sarimatest.c)

target_link_libraries(sarimatest ctsalib m)

add_executable(misctests misctests.c)

target_link_libraries(misctests ctsalib m)


set_target_properties(acftest arimatest artest artest2 pacftest sarimatest misctests
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
)
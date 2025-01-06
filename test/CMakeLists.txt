# SPDX-License-Identifier: BSD-3-Clause
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

# TODO add the header files via a `target_include_directories`
add_executable(sarimatest2 
			   sarimatest2.c 
			   auxlib.c 
			   auxlib.h 
			   csv_in.c 
			   csv_in.h 
			   )

target_link_libraries(sarimatest2 ctsalib m)

# TODO add the header files via a `target_include_directories`
add_executable(sarimatest3 
			   sarimatest3.c 
			   auxlib.c 
			   auxlib.h 
			   csv_in.c 
			   csv_in.h 
			   )

target_link_libraries(sarimatest3 ctsalib m)

add_executable(misctests misctests.c)

target_link_libraries(misctests ctsalib m)

add_executable(dwttests dwttests.c)

target_link_libraries(dwttests ctsalib m)

add_executable(sarimaxtest1 sarimaxtest1.c)

target_link_libraries(sarimaxtest1 ctsalib m)

add_executable(sarimaxtest2 sarimaxtest2.c)

target_link_libraries(sarimaxtest2 ctsalib m)

add_executable(sarimaxtest3 sarimaxtest3.c)

target_link_libraries(sarimaxtest3 ctsalib m)

add_executable(autoarimatest1 autoarimatest1.c)

target_link_libraries(autoarimatest1 ctsalib m)

add_executable(autoarimatest2 autoarimatest2.c)

target_link_libraries(autoarimatest2 ctsalib m)

add_executable(autoarimatest3 autoarimatest3.c)

target_link_libraries(autoarimatest3 ctsalib m)


set_target_properties(acftest arimatest artest artest2 pacftest sarimatest misctests dwttests sarimaxtest1 sarimaxtest2 sarimaxtest3 autoarimatest1 autoarimatest2 autoarimatest3
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}"
)

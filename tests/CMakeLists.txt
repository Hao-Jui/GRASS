# Unit tests
set(UNIT_TESTS test_ad test_spline test_spectral test_brent test_eos)
foreach(test ${UNIT_TESTS})
  add_executable(${test} unit/${test}.f90 unit/test_utils.f90)
  target_link_libraries(${test} grass_lib)
  add_test(NAME ${test} COMMAND ${test}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()

# Integration tests
set(INTEG_TESTS test_gr_uniform test_gr_constj test_gr_uryu
  test_st_uniform test_st_uniform_r07)
foreach(test ${INTEG_TESTS})
  add_executable(${test} integration/${test}.f90)
  target_link_libraries(${test} grass_lib)
  add_test(NAME ${test} COMMAND ${test}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  set_tests_properties(${test} PROPERTIES TIMEOUT 120)
endforeach()

# Ensure output directories exist for integration/restart tests
add_test(NAME create_output_dirs
  COMMAND ${CMAKE_COMMAND} -E make_directory
    ${CMAKE_SOURCE_DIR}/Cont ${CMAKE_SOURCE_DIR}/Res)
set_tests_properties(create_output_dirs PROPERTIES FIXTURES_SETUP OutputDirs)
set_tests_properties(${INTEG_TESTS} PROPERTIES FIXTURES_REQUIRED OutputDirs)

# Restart test
add_executable(test_restart integration/test_restart.f90)
target_link_libraries(test_restart grass_lib)
add_test(NAME test_restart COMMAND test_restart
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(test_restart PROPERTIES TIMEOUT 120 FIXTURES_REQUIRED OutputDirs)

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

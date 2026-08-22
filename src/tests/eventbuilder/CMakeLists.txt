find_package(GTest QUIET)
if(NOT GTest_FOUND)
  message(STATUS "GTest not found — skipping eventbuilder tests")
  return()
endif()

get_filename_component(TEST_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

add_executable(${TEST_NAME} EventBuilderTest.cpp)

target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(TARGET onnxruntime::onnxruntime)
  target_compile_definitions(${TEST_NAME} PRIVATE HAS_ONNXRUNTIME)
  target_link_libraries(${TEST_NAME}
    PRIVATE GTest::gtest GTest::gtest_main onnxruntime::onnxruntime)
else()
  target_link_libraries(${TEST_NAME}
    PRIVATE GTest::gtest GTest::gtest_main)
endif()

install(TARGETS ${TEST_NAME} DESTINATION bin)

add_test(NAME t_${TEST_NAME} COMMAND $<TARGET_FILE:${TEST_NAME}>)

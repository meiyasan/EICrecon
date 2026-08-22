set(PLUGIN_NAME "onnx_service")

plugin_add(${PLUGIN_NAME} WITH_SHARED_LIBRARY WITHOUT_PLUGIN)
plugin_glob_all(${PLUGIN_NAME})
plugin_add_onnxruntime(${PLUGIN_NAME})

# nanoflann v1.6.2 vendored as a single header (nanoflann.hpp in this directory).
target_compile_definitions(${PLUGIN_NAME}_library PRIVATE USE_NANOFLANN)
target_include_directories(${PLUGIN_NAME}_library PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

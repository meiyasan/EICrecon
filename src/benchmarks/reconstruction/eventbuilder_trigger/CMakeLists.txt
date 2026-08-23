# Automatically set plugin name the same as the directory name
get_filename_component(PLUGIN_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

plugin_add(${PLUGIN_NAME})

plugin_add_event_model(${PLUGIN_NAME})
plugin_add_cern_root(${PLUGIN_NAME})

plugin_glob_all(${PLUGIN_NAME})

plugin_link_libraries(${PLUGIN_NAME} log_library ROOT::Hist ROOT::Gpad ROOT::Graf)

cmake_minimum_required(VERSION 3.16)

if(JANA_VERSION VERSION_LESS "2.4.3")

  message(STATUS "Skipping `eventbuilder` plugin because JANA_VERSION < 2.4.3")

else()

  get_filename_component(PLUGIN_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

  # Function creates ${PLUGIN_NAME}_plugin and ${PLUGIN_NAME}_library targets
  # Setting default includes, libraries and installation paths
  plugin_add(${PLUGIN_NAME} PLUGIN_USE_CC_ONLY)

  # The macro grabs sources as *.cc *.cpp *.c and headers as *.h *.hh *.hpp Then
  # correctly sets sources for ${_name}_plugin and ${_name}_library targets Adds
  # headers to the correct installation directory
  plugin_glob_all(${PLUGIN_NAME})

  # plugin_glob_all only globs this directory; pull in the per-detector
  # digitization registration files that live in the digi/ subdirectory.
  file(GLOB DIGITS_SRC CONFIGURE_DEPENDS ${CMAKE_CURRENT_LIST_DIR}/plugins/*.cc)
  target_sources(${PLUGIN_NAME}_plugin PRIVATE ${DIGITS_SRC})
  if(TARGET ${PLUGIN_NAME}_library)
    target_sources(${PLUGIN_NAME}_library PRIVATE ${DIGITS_SRC})
  endif()

  # Find dependencies plugin_add_dd4hep(${PLUGIN_NAME}#)
  # plugin_add_acts(${PLUGIN_NAME}) plugin_add_cern_root(${PLUGIN_NAME})
  plugin_add_event_model(${PLUGIN_NAME})
  plugin_add_onnxruntime(${PLUGIN_NAME})

  # Add include directories (works same as target_include_directories)
  # plugin_include_directories(${PLUGIN_NAME} SYSTEM PUBLIC ... )

  # Add libraries (same as target_include_directories but for both plugin and
  # library)
  plugin_link_libraries(${PLUGIN_NAME} algorithms_digi_library
                        algorithms_tracking_library
                        onnx_service_library)
  #
  # Add include directories (works same as target_include_directories)
  # plugin_include_directories(${PLUGIN_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR}
  # SYSTEM PUBLIC ${podio_INCLUDE_DIR} ${EDM4HEP_INCLUDE_DIR}
  # ${DD4hep_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${EDM4EIC_INCLUDE_DIRS})
  #
  # Add libraries (works same as target_include_directories)
  # plugin_link_libraries(${PLUGIN_NAME} ${ROOT_LIBRARIES}
  # algorithms_tracking_library EDM4HEP::edm4hep EDM4EIC::edm4eic spdlog::spdlog
  # Boost::boost ${ROOT_EG_LIBRARY})

endif()

# ── test-events ────────────────────────────────────────────────────────────────
# Run eicrecon+eventbuilder on 1 edm4hep event per SRO class.
# Writes <class>.eicrecon.root to vendor/eicrecon/share/samples/gold-new/.
#
# Usage inside the container:
#   make test-events
#
# The target searches vendor/sro/datasets/edm4hep/{vanilla,gold-new,gold-old}/<class>/
# for input files and skips any class for which no edm4hep input is present.
# Run `make dataset DATASET=vanilla` in vendor/sro first to populate inputs.
#
# Override the eicrecon executable:
#   cmake --build . --target test-events -- EICRECON_BIN=/path/to/eicrecon

set(_EB_SAMPLES_DIR "${CMAKE_SOURCE_DIR}/vendor/eicrecon/share/samples/gold-new")
set(_EB_EDM4HEP_ROOT "${CMAKE_SOURCE_DIR}/vendor/sro/datasets/edm4hep")
set(_EB_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/scripts/test_events.sh")

if(NOT DEFINED EICRECON_BIN OR EICRECON_BIN STREQUAL "")
  find_program(_EB_EICRECON_BIN eicrecon)
  if(_EB_EICRECON_BIN)
    set(EICRECON_BIN "${_EB_EICRECON_BIN}")
  else()
    set(EICRECON_BIN "eicrecon")  # fall back to PATH at run time
  endif()
endif()

add_custom_target(test-events
  COMMAND bash "${_EB_SCRIPT}"
               "${_EB_SAMPLES_DIR}"
               "${_EB_EDM4HEP_ROOT}"
               "${EICRECON_BIN}"
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  COMMENT "Running eicrecon+eventbuilder on 1 event per SRO class → share/samples/gold-new/"
  VERBATIM
)

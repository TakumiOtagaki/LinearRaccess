# Re-evaluate at every build, including after commits without source edits.
# configure_file only changes the header when its content changes.
set(SOURCE_REVISION "unknown")
find_program(GIT_EXECUTABLE git)
if(GIT_EXECUTABLE AND EXISTS "${SOURCE_DIR}/.git")
  execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${SOURCE_DIR}" rev-parse HEAD
    OUTPUT_VARIABLE GIT_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE GIT_RESULT ERROR_QUIET)
  execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${SOURCE_DIR}" status --porcelain --untracked-files=normal
    OUTPUT_VARIABLE GIT_STATUS OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE STATUS_RESULT ERROR_QUIET)
  if(GIT_RESULT EQUAL 0 AND STATUS_RESULT EQUAL 0 AND GIT_REVISION MATCHES "^[0-9a-f]+$")
    set(SOURCE_REVISION "${GIT_REVISION}")
    if(NOT GIT_STATUS STREQUAL "")
      string(APPEND SOURCE_REVISION "-dirty")
    endif()
  endif()
endif()
file(MAKE_DIRECTORY "${OUTPUT_DIR}")
configure_file("${CMAKE_CURRENT_LIST_DIR}/build_revision.hpp.in"
               "${OUTPUT_DIR}/build_revision.hpp" @ONLY)

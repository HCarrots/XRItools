cmake_minimum_required(VERSION 3.24)
if(POLICY CMP0207)
    cmake_policy(SET CMP0207 NEW)
endif()
if(NOT DEFINED APP OR NOT DEFINED OUTPUT OR NOT DEFINED SEARCH_PATHS)
    message(FATAL_ERROR "APP, OUTPUT and SEARCH_PATHS must be supplied")
endif()
file(GLOB_RECURSE plugin_dlls "${OUTPUT}/platforms/*.dll" "${OUTPUT}/styles/*.dll"
    "${OUTPUT}/imageformats/*.dll" "${OUTPUT}/iconengines/*.dll")
file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES "${APP}"
    LIBRARIES ${plugin_dlls}
    DIRECTORIES ${SEARCH_PATHS}
    PRE_EXCLUDE_REGEXES "api-ms-.*" "ext-ms-.*"
    POST_EXCLUDE_REGEXES ".*[/\\\\][Ss][Yy][Ss][Tt][Ee][Mm]32[/\\\\].*"
    CONFLICTING_DEPENDENCIES_PREFIX conflicts
    RESOLVED_DEPENDENCIES_VAR resolved
    UNRESOLVED_DEPENDENCIES_VAR unresolved)
foreach(filename IN LISTS conflicts_FILENAMES)
    set(preferred "")
    foreach(directory IN LISTS SEARCH_PATHS)
        if(EXISTS "${directory}/${filename}")
            set(preferred "${directory}/${filename}")
            break()
        endif()
    endforeach()
    if(NOT preferred)
        message(FATAL_ERROR "Conflicting dependency without an SDK candidate: ${filename}")
    endif()
    list(APPEND resolved "${preferred}")
endforeach()
foreach(dependency IN LISTS resolved)
    string(TOLOWER "${dependency}" lower_path)
    if(NOT lower_path MATCHES "^[a-z]:/windows/")
        file(COPY "${dependency}" DESTINATION "${OUTPUT}")
    endif()
endforeach()
if(unresolved)
    message(FATAL_ERROR "Unresolved runtime dependencies: ${unresolved}")
endif()
message(STATUS "Runtime dependencies copied to ${OUTPUT}")

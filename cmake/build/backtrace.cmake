

set (backtrace_src
        ../../third-party/backtrace/backtrace.c
        )
set_source_files_properties(
            backtrace_src
            PROPERTIES
            LANGUAGE C
            )
MESSAGE(${COMPILE_FLAGS})
set(OLD_FLAGS ${COMPILE_FLAGS})

set(COMPILE_FLAGS "-O2 -shared -Wall -lbfd -liberty -limagehlp -lz -lintl")

MESSAGE(${COMPILE_FLAGS})
add_library(backtrace ${backtrace_src})
target_link_libraries(backtrace_src)

set(COMPILE_FLAGS OLD_FLAGS)

MESSAGE(COMPILE_FLAGS)

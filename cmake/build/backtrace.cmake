

set (backtrace_src
        ../../third-party/backtrace/backtrace.c
        )
set_source_files_properties(
            backtrace_src
            PROPERTIES
            LANGUAGE C
            )

set(OLD_FLAGS ${COMPILE_FLAGS})

set(COMPILE_FLAGS "-O2 -shared -o -Wall -lbfd -liberty -limagehlp -lz -lintl")


add_library(backtrace ${backtrace_src})

target_link_libraries(backtrace)

set(COMPILE_FLAGS ${OLD_FLAGS})


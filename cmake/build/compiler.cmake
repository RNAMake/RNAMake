if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(COMPILER clang)
  # using Clang
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(COMPILER gcc)
  # using GCC
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(COMPILER AppleClang)
    #endif()
    #
    ## what do we want to do here? 
    #if( ${CMAKE_CXX_COMPILER} MATCHES "[g][+][+].*" )
    #elseif(${CMAKE_CXX_COMPILER} MATCHES ".*clang[+]{2}.*")
else()
    message(FATAL_ERROR
        "${BoldRed}Specified C++ compiler is \"${CMAKE_CXX_COMPILER}\". "
        "On this machine, the variable \"CMAKE_CXX_COMPILER\" is ambiguous. "
        "RNAMake's build currently supports clang and gcc.${ColourReset}\n"
        "${BoldGreen}For successful build, try the command:\n"
        "\t$ cmake [-G Ninja] -DCMAKE_CXX_COMPILER=[g++|clang].${ColourReset}\n\n" )
endif()

MESSAGE(
"${BoldGreen}>> CMAKE identifies C++ compiler as '${CMAKE_CXX_COMPILER}',\n\t\t==> interpreting this as '${COMPILER}'\n" 
">> To change, set CXX and CC environment variables (or pass -DCMAKE_CXX_COMPILER) and do a clean rebuild.\n"
">>  current settings: CXX='$ENV{CXX}' CC='$ENV{CC}'${ColourReset}\n" 
)

if( ${COMPILER} STREQUAL clang  OR ${COMPILER} STREQUAL AppleClang  )
	set(compile_flags 
		-Wno-c++11-extensions;
        -std=gnu++1y;
	)


endif()

if( ${COMPILER} STREQUAL gcc )
	set(compile_flags
	-std=c++14;
    	-w;
    	-O3;
    	-pedantic;
	)
    

endif()

foreach( flag ${compile_flags} )
    set( COMPILE_FLAGS "${COMPILE_FLAGS} ${flag}" )
endforeach()

foreach( flag ${COMPILE_FLAGS} )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" )
endforeach()


if(CMAKE_SYSTEM_NAME STREQUAL "Windows")  
    include(backtrace.cmake) 

elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread " )
    set( COMPILE_FLAGS "${COMPILE_FLAGS} -pthread " )
    set( CMAKE_CXX_LINKER_FLAGS "${CMAKE_CXX_LINKER_FLAGS} -ldl " )
    set( CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-as-needed -ldl")
    set( CMAKE_EXE_LINKER_FLAGS " -lstdc++ -Wl,--no-as-needed,--no-export-dynamic ")

elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")

    SET(CMAKE_C_ARCHIVE_CREATE   "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
    SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
    SET(CMAKE_C_ARCHIVE_FINISH   "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
    SET(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")

endif()


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
		-Wno-c++11-extensions
        -std=gnu++1y
        -g
        -pedantic
        #    -O4 
	)
endif()

if( ${COMPILER} STREQUAL gcc )
	set(compile_flags
	-std=c++14
    -w
    #-O3
    -pedantic
    -Wall
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
endif()

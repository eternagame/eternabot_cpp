if( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" )
    set(COMPILER clang)
else()
    set(COMPILER gcc)
endif()

MESSAGE( ">> CMAKE identifies C++ compiler as '${CMAKE_CXX_COMPILER_ID}', interpreting this as '${COMPILER}'" )
MESSAGE( ">> To change, set CXX and CC environment variables (or pass -DCMAKE_CXX_COMPILER) and do a clean rebuild.")
MESSAGE( ">>  current settings: CXX='$ENV{CXX}' CC='$ENV{CC}'" )

if( ${COMPILER} STREQUAL clang )
	set(compile_flags 
		-Wno-c++11-extensions
        -std=gnu++1y
        -g
	)
endif()

if( ${COMPILER} STREQUAL gcc )
	set(compile_flags
	-std=c++14
    -w
	)
endif()

foreach( flag ${compile_flags} )
    set( COMPILE_FLAGS "${COMPILE_FLAGS} ${flag}" )
endforeach()

foreach( flag ${COMPILE_FLAGS} )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" )
endforeach()

MESSAGE(“${FLAGS}”)



cmake_minimum_required(VERSION 2.8.12)
project(rnamake_new)

set(CMAKE_BUILD_TYPE Release)
include(/Users/josephyesselman/projects/eternabot_cpp/cmake/build/compiler.cmake)

# Include path for Python header files
# Include paths for RNAMake src
include_directories(/Users/josephyesselman/projects/eternabot_cpp/src/)
include_directories(/Users/josephyesselman/projects/eternabot_cpp/src/plog/)
include_directories(/Users/josephyesselman/projects/eternabot_cpp/unittests/)

include_directories(/Users/josephyesselman/projects/eternabot_cpp/apps/)

###########################################################################
# base
###########################################################################

set(base_files
	/Users/josephyesselman/projects/eternabot_cpp/src/base/file_io.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/base/string.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/base/application.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/base/cl_option.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/base/command_line_parser.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/base/backtrace.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/base/log.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/base/option.cc
)

add_library(base_lib ${base_files})
target_link_libraries(base_lib)

###########################################################################
# math
###########################################################################

set(math_files
	/Users/josephyesselman/projects/eternabot_cpp/src/math/math.cc
)

add_library(math_lib ${math_files})
target_link_libraries(math_lib base_lib )

###########################################################################
# data_structure
###########################################################################

set(data_structure_files
	/Users/josephyesselman/projects/eternabot_cpp/src/data_structure/graph/graph.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/data_structure/graph/graph_node.cc
)

add_library(data_structure_lib ${data_structure_files})
target_link_libraries(data_structure_lib base_lib )

###########################################################################
# util
###########################################################################

set(util_files
	/Users/josephyesselman/projects/eternabot_cpp/src/util/motif_type.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/util/uuid.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/util/monte_carlo.cpp
)

add_library(util_lib ${util_files})
target_link_libraries(util_lib math_lib )

###########################################################################
# vienna
###########################################################################

set(vienna_files
	/Users/josephyesselman/projects/eternabot_cpp/src/vienna/energy_par.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/vienna/vienna.cc
)

add_library(vienna_lib ${vienna_files})
target_link_libraries(vienna_lib base_lib )

###########################################################################
# secondary_structure
###########################################################################

set(secondary_structure_files
	/Users/josephyesselman/projects/eternabot_cpp/src/secondary_structure/sequence_tools.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/secondary_structure/util.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/secondary_structure/structure.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/secondary_structure/residue.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/secondary_structure/motif.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/secondary_structure/basepair.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/secondary_structure/rna_structure.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/secondary_structure/chain.cc
	/Users/josephyesselman/projects/eternabot_cpp/src/secondary_structure/pose.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/secondary_structure/secondary_structure_parser.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/secondary_structure/sequence_constraint.cc
)

add_library(secondary_structure_lib ${secondary_structure_files})
target_link_libraries(secondary_structure_lib util_lib )

###########################################################################
# eternabot
###########################################################################

set(eternabot_files
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/sequence_designer.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/scorer.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/strategy.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/feature_generator.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/strategy/modified_num_of_yellow.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/strategy/num_of_yellow.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/strategy/modified_clear_plot.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/strategy/direction_of_gc.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/strategy/berex_test.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/strategy/modified_direction_of_gc.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/strategy/modified_berex_test.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/strategy/clear_plot.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/strategy/modified_a_basic_test.cpp
	/Users/josephyesselman/projects/eternabot_cpp/src/eternabot/strategy/a_basic_test.cpp
)

add_library(eternabot_lib ${eternabot_files})
target_link_libraries(eternabot_lib vienna_lib  secondary_structure_lib )

###########################################################################
# all
###########################################################################

add_library(all_lib /Users/josephyesselman/projects/eternabot_cpp/src/main.cpp)
target_link_libraries(all_lib eternabot_lib )

add_executable(eternabot /Users/josephyesselman/projects/eternabot_cpp/apps/eternabot/eternabot.cpp )
target_link_libraries(eternabot all_lib)


add_executable(score_rnas /Users/josephyesselman/projects/eternabot_cpp/apps/score_rnas/score_rnas.cc )
target_link_libraries(score_rnas all_lib)



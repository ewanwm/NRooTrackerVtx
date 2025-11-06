
project(nRooTrackerConverter)

find_package(ROOT)

message("${CMAKE_CURRENT_SOURCE_DIR}")

# Generate ROOT dictionaries for NRooTrackerVtx classes
ROOT_GENERATE_DICTIONARY(
        TJNuBeamFlux_dict
        ${CMAKE_CURRENT_SOURCE_DIR}/src/TJNuBeamFlux.hxx              
        LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/src/TJNuBeamFlux_LinkDef.hxx
)


ROOT_GENERATE_DICTIONARY(
	TNRooTrackerVtx_dict 
	${CMAKE_CURRENT_SOURCE_DIR}/src/TNRooTrackerVtx.hxx 
	LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/src/TNRooTrackerVtx_LinkDef.hxx
)

# Create NRooTrackerVtx library
add_library(
	NRooTrackerVtx 
	SHARED 
	TJNuBeamFlux_dict.cxx 
	TNRooTrackerVtx_dict.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/src/TNRooTrackerVtx.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/src/TJNuBeamFlux.cxx
)

target_link_libraries(
	NRooTrackerVtx
	PUBLIC
	ROOT::Core
	ROOT::RIO
	ROOT::Tree
)

target_include_directories(
        NRooTrackerVtx
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Make executables
add_executable(
	NRooTrackerConverter
	${CMAKE_CURRENT_SOURCE_DIR}/app/NRooTrackerConverter.cxx
)

target_link_libraries(NRooTrackerConverter NRooTrackerVtx ROOT::Core)


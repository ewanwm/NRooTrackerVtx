
project(nRooTrackerConverter)

find_package(ROOT)

message("${CMAKE_CURRENT_SOURCE_DIR}")

# Generate ROOT dictionaries for TJNuBeamFlux class
add_library(
	TJNuBeamFlux_dict
	SHARED
)

target_link_libraries(
	TJNuBeamFlux_dict
	PUBLIC
	${ROOT_LIBRARIES}
	ROOT::Core
	ROOT::RIO
	ROOT::Tree
)

ROOT_GENERATE_DICTIONARY(
	TJNuBeamFlux_dict
	${CMAKE_CURRENT_SOURCE_DIR}/src/TJNuBeamFlux.hxx
	LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/src/TJNuBeamFlux_LinkDef.hxx
	MODULE TJNuBeamFlux_dict
)

target_sources(
	TJNuBeamFlux_dict
	PUBLIC
	${CMAKE_CURRENT_BINARY_DIR}/TJNuBeamFlux_dict.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/src/TJNuBeamFlux.cxx
)


# Generate ROOT dictionaries for TNRooTrackerVtx class
add_library(
	TNRooTrackerVtx_dict
	SHARED
)

target_link_libraries(
	TNRooTrackerVtx_dict
	PUBLIC
	${ROOT_LIBRARIES}
	ROOT::Core
	ROOT::RIO
	ROOT::Tree
)

ROOT_GENERATE_DICTIONARY(
	TNRooTrackerVtx_dict 
	${CMAKE_CURRENT_SOURCE_DIR}/src/TNRooTrackerVtx.hxx 
	LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/src/TNRooTrackerVtx_LinkDef.hxx
	MODULE TNRooTrackerVtx_dict
)

target_sources(
	TNRooTrackerVtx_dict
	PUBLIC
	${CMAKE_CURRENT_BINARY_DIR}/TJNuBeamFlux_dict.cxx
	${CMAKE_CURRENT_BINARY_DIR}/TNRooTrackerVtx_dict.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/src/TNRooTrackerVtx.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/src/TJNuBeamFlux.cxx
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
	TNRooTrackerVtx_dict
	TJNuBeamFlux_dict
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


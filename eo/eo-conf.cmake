
# Current version
SET(PROJECT_VERSION_MAJOR 1)
SET(PROJECT_VERSION_MINOR 3)
SET(PROJECT_VERSION_PATCH 2)
SET(PROJECT_VERSION_MISC "-edge")

# ADD_DEFINITIONS(-DDEPRECATED_MESSAGES) # disable warning deprecated function messages
# If you plan to use OpenMP, put the following boolean to true :
SET(WITH_OMP FALSE CACHE BOOL "Use OpenMP ?" FORCE)

# If you plan to use MPI, precise here where are the static libraries from
# openmpi and boost::mpi.

SET(WITH_MPI TRUE CACHE BOOL "Use mpi ?" FORCE)
SET(WITH_BOOST TRUE CACHE BOOL "Use boost ?" FORCE)
SET(MPI_DIR "/usr/local" CACHE PATH "OpenMPI directory" FORCE)

SET(CXX11 TRUE CACHE BOOL "Use c++11 ?" FORCE)

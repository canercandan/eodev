# File: FindEO.cmake
# CMAKE commands to actually use the EO library
# Version: 0.0.1
#
# The following variables are filled out:
# - EO_INCLUDE_DIRS
# - EO_LIBRARY_DIRS
# - EO_LIBRARIES
# - EO_FOUND
#
# Here are the components:
# - PyEO
# - es
# - ga
# - cma
# - mpi
# - serial
#
# You can use FIND_PACKAGE( EO COMPONENTS ... ) to enable one or several components.
#

# Default enabled components
SET(EO_LIBRARIES_TO_FIND eo eoutils)

# Use FIND_PACKAGE( EO COMPONENTS ... ) to enable modules
IF(EO_FIND_COMPONENTS)
  FOREACH(component ${EO_FIND_COMPONENTS})
    STRING(TOUPPER ${component} _COMPONENT)
    SET(EO_USE_${_COMPONENT} 1)
  ENDFOREACH(component)

  # To make sure we don't use PyEO, ES, GA, CMA, MPI, SERIAL when there arent in COMPONENTS
  IF(NOT EO_USE_PYEO)
    SET(EO_DONT_USE_PYEO 1)
  ELSE(NOT EO_USE_PYEO)
    SET(EO_LIBRARIES_TO_FIND ${EO_LIBRARIES_TO_FIND} PyEO)
  ENDIF(NOT EO_USE_PYEO)

  IF(NOT EO_USE_ES)
    SET(EO_DONT_USE_ES 1)
  ELSE(NOT EO_USE_ES)
    SET(EO_LIBRARIES_TO_FIND ${EO_LIBRARIES_TO_FIND} es)
  ENDIF(NOT EO_USE_ES)

  IF(NOT EO_USE_GA)
    SET(EO_DONT_USE_GA 1)
  ELSE(NOT EO_USE_GA)
    SET(EO_LIBRARIES_TO_FIND ${EO_LIBRARIES_TO_FIND} ga)
  ENDIF(NOT EO_USE_GA)

  IF(NOT EO_USE_CMA)
    SET(EO_DONT_USE_CMA 1)
  ELSE(NOT EO_USE_CMA)
    SET(EO_LIBRARIES_TO_FIND ${EO_LIBRARIES_TO_FIND} cma)
  ENDIF(NOT EO_USE_CMA)

  IF(NOT EO_USE_MPI)
    SET(EO_DONT_USE_MPI 1)
  ELSE(NOT EO_USE_MPI)
    SET(EO_LIBRARIES_TO_FIND ${EO_LIBRARIES_TO_FIND} eompi)
  ENDIF(NOT EO_USE_MPI)

  IF(NOT EO_USE_SERIAL)
    SET(EO_DONT_USE_SERIAL 1)
  ELSE(NOT EO_USE_SERIAL)
    SET(EO_LIBRARIES_TO_FIND ${EO_LIBRARIES_TO_FIND} eoserial)
  ENDIF(NOT EO_USE_SERIAL)
ENDIF(EO_FIND_COMPONENTS)

IF(NOT EO_INCLUDE_DIRS)
  FIND_PATH(
    EO_INCLUDE_DIRS
    EO.h
    PATHS
    /usr/include/eo
    /usr/local/include/eo
    )
ENDIF(NOT EO_INCLUDE_DIRS)

IF(NOT EO_LIBRARY_DIRS)
  FIND_PATH(
    EO_LIBRARY_DIRS
    libeo.a
    PATHS
    /usr/lib
    /usr/local/lib
    )
ENDIF(NOT EO_LIBRARY_DIRS)

IF(NOT EO_LIBRARIES)
  SET(EO_LIBRARIES)
  FOREACH(component ${EO_LIBRARIES_TO_FIND})
    FIND_LIBRARY(
      EO_${component}_LIBRARY
      NAMES ${component}
      PATHS
      /usr/lib
      /usr/local/lib
      )
    IF(EO_${component}_LIBRARY)
      SET(EO_LIBRARIES ${EO_LIBRARIES} ${EO_${component}_LIBRARY})
    ELSE(EO_${component}_LIBRARY)
      MESSAGE(FATAL_ERROR "${component} component not found.")
    ENDIF(EO_${component}_LIBRARY)
  ENDFOREACH(component)
ENDIF(NOT EO_LIBRARIES)

# handle the QUIETLY and REQUIRED arguments and set EO_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(EO
  REQUIRED_VARS EO_INCLUDE_DIRS EO_LIBRARY_DIRS EO_LIBRARIES
  )

MARK_AS_ADVANCED(EO_INCLUDE_DIRS)
MARK_AS_ADVANCED(EO_LIBRARY_DIRS)
MARK_AS_ADVANCED(EO_LIBRARIES)

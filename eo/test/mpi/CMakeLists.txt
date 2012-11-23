###############################################################################
##
## CMakeLists file for eo/test/mpi
##
###############################################################################

######################################################################################
### 1) Include the sources
######################################################################################

MESSAGE("EO SOURCE DIR: ${EO_SOURCE_DIR}")
MESSAGE("OMPI: ${MPI_DIR}")

INCLUDE_DIRECTORIES(${MPI_DIR}/include)
INCLUDE_DIRECTORIES(${EO_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(/Dev/paradiseo/problems/src)

######################################################################################
### 2) Specify where CMake can find the libraries
######################################################################################

LINK_DIRECTORIES(${EO_BINARY_DIR}/lib)
LINK_DIRECTORIES(${MPI_DIR}/lib)

######################################################################################
### 3) Define your targets and link the librairies
######################################################################################

SET (TEST_LIST
    t-mpi-parallelApply
    t-mpi-wrapper
    t-mpi-multipleRoles
    t-mpi-eval
    t-mpi-multistart
    t-mpi-distrib-exp
    # t-mpi-justatest
    # t-mpi-eopop
    # t-mpi-sendRecvPop
    # t-mpi-sendRecvPop-with-thread
    # t-mpi-sendThreadRecvThreadPop
    # t-mpi-allSendRecvPop-with-thread
    # t-mpi-allSendRecvPop-with-thread-2
    # t-mpi-allSendRecvPop-with-thread-3
    # t-mpi-allSendRecvPop-with-thread-4
    # t-mpi-allSendRecvPop-with-thread-4-serial
    # t-mpi-isendIrecvPop
    # t-mpi-island-model
    #t-mpi-thread
    # t-mpi-im-with-thread
    # t-mpi-im-with-thread-inverted
    # t-thread
    # t-thread-notify
    # t-thread-exit
    # t-mpi-nb
    # t-mpi-graph
    # t-mpi-graph-sendrecv
    # t-mpi-graph-isendirecv
    # t-mpi-graph-isendirecv_eompi
    # t-mpi-graph-isendirecv_eompi_pop_alone
    # t-mpi-graph-isendirecv_eompi_pop_alone2
    # t-mpi-graph-isendirecv_eompi_pop_alone_with_thread
    # t-mpi-graph-isendirecv_eompi_pop_alone_with_thread2
    # t-mpi-graph-isendirecv_eompi_pop_alone_with_thread3
    # t-mpi-unpack
    )

FOREACH (test ${TEST_LIST})
  SET ("T_${test}_SOURCES" "${test}.cpp")
ENDFOREACH (test)

SET(CMAKE_CXX_COMPILER "mpicxx")
ADD_DEFINITIONS(-DWITH_MPI)

# IF(ENABLE_CMAKE_TESTING)
    # ADD_LIBRARY(mpi-thread thread.cpp)
    FOREACH (test ${TEST_LIST})
      ADD_EXECUTABLE(${test} ${T_${test}_SOURCES})
      ADD_TEST(${test} ${test})
      TARGET_LINK_LIBRARIES(${test} eompi eoserial ga eoutils eo)
      INSTALL(TARGETS ${test} RUNTIME DESTINATION share/eo/test COMPONENT test)
    ENDFOREACH (test)
# ENDIF()

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

SET (BOOST_TEST_LIST
    t-boost-mpi-graph-isendirecv_eompi_pop_alone2
    t-boost-mpi-graph-isendirecv_eompi_pop_alone2_mt
    )

FOREACH (test ${BOOST_TEST_LIST})
  SET ("BOOST_T_${test}_SOURCES" "${test}.cpp")
ENDFOREACH (test)

# IF(ENABLE_CMAKE_TESTING)
    FOREACH (test ${BOOST_TEST_LIST})
      ADD_EXECUTABLE(${test} ${BOOST_T_${test}_SOURCES})
      ADD_TEST(${test} ${test})
      TARGET_LINK_LIBRARIES(${test} eoutils ${Boost_LIBRARIES} eompi eoserial ga eoutils eo)
      INSTALL(TARGETS ${test} RUNTIME DESTINATION share/eo/test COMPONENT test)
    ENDFOREACH (test)
# ENDIF()

SET(RESOURCES
  gatherMonitor.py
  )

FOREACH(file ${RESOURCES})
  EXECUTE_PROCESS(
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/${file}
    ${CMAKE_CURRENT_BINARY_DIR}/${file}
    )
ENDFOREACH(file)

######################################################################################

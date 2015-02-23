set(SOURCES
  main.cpp
  CartesianProduct.cpp
  ExpMapMatrix.cpp
  Manifold.cpp
  Point.cpp
  RealSpace.cpp
  ReusableTemporaryMap.cpp
  utils.cpp
  )
set(HEADERS
  ../include/manifolds/CartesianProduct.h
  ../include/manifolds/defs.h
  ../include/manifolds/ExpMapMatrix.h
  ../include/manifolds/Manifold.h
  ../include/manifolds/pgs_assert.h
  ../include/manifolds/Point.h
  ../include/manifolds/RealSpace.h
  ../include/manifolds/ReusableTemporaryMap.h
  ../include/manifolds/SO3.h
  ../include/manifolds/utils.h
  ../include/manifolds/view.h
  ../include/manifolds/RefCounter.h
  ../include/manifolds/ValidManifold.h
  ../include/manifolds/manifolds_api.h
)

# set(LIBRARIES
#   manifolds
#   manifoldsTest
# )
# foreach(lib ${LIBRARIES})
#   add_library(${lib} SHARED ${SOURCES} ${HEADERS})
#   set_target_properties(${lib} PROPERTIES
#     VERSION ${PROJECT_VERSION}
#     SOVERSION 0.1
#   )
#   target_compile_definitions(${lib} PRIVATE "-DMANIFOLDS_EXPORT")
# endforeach()

# manifolds library
add_library(manifolds SHARED ${SOURCES} ${HEADERS})
set_target_properties(manifolds PROPERTIES
	                  VERSION ${PROJECT_VERSION}
                      SOVERSION 0.1)
target_compile_definitions(manifolds PRIVATE "-DMANIFOLDS_EXPORT")

# manifoldsTest library
add_library(manifoldsTest SHARED ${SOURCES} ${HEADERS})
set_target_properties(manifoldsTest PROPERTIES
	                  VERSION ${PROJECT_VERSION}
                      SOVERSION 0.1)
#target_compile_definitions(manifoldsTest PRIVATE -DMANIFOLDS_EXPORT -DPGS_ASSERT_THROW -DEIGEN_RUNTIME_NO_MALLOC)
target_compile_definitions(manifoldsTest PRIVATE -DMANIFOLDS_EXPORT -DPGS_ASSERT_THROW)

install(TARGETS manifolds
  RUNTIME DESTINATION ${CMAKE_INSTALL_FULL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
)
if(WIN32)
install(TARGETS manifoldsTest
  RUNTIME DESTINATION ${CMAKE_INSTALL_FULL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
)
endif()

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/manifolds)

add_executable(ManifoldsExe main.cpp)
target_link_libraries(ManifoldsExe manifolds)

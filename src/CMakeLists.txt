## Copyright (c) 2015 CNRS
## Authors: Stanislas Brossette, Adrien Escande, Pierre Gergondet 

## This file is part of manifolds
## manifolds is free software: you can redistribute it
## and/or modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation, either version
## 3 of the License, or (at your option) any later version.

## manifolds is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied warranty
## of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Lesser Public License for more details.  You should have
## received a copy of the GNU Lesser General Public License along with
## manifolds. If not, see
## <http://www.gnu.org/licenses/>.

set(SOURCES
  CartesianProduct.cpp
  CartesianPower.cpp
  ExpMapMatrix.cpp
  ExpMapQuaternion.cpp
  Manifold.cpp
  Point.cpp
  RealSpace.cpp
  ReusableTemporaryMap.cpp
  S2.cpp
  utils.cpp
  )
set(HEADERS
  ../include/manifolds/CartesianProduct.h
  ../include/manifolds/CartesianPower.h
  ../include/manifolds/defs.h
  ../include/manifolds/ExpMapMatrix.h
  ../include/manifolds/ExpMapQuaternion.h
  ../include/manifolds/Manifold.h
  ../include/manifolds/mnf_assert.h
  ../include/manifolds/Point.h
  ../include/manifolds/RealSpace.h
  ../include/manifolds/ReusableTemporaryMap.h
  ../include/manifolds/SO3.h
  ../include/manifolds/S2.h
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
set_target_properties(manifolds PROPERTIES
                    COMPILE_DEFINITIONS "PRIVATE -DMANIFOLDS_EXPORT")

# manifoldsTest library
add_library(manifoldsTest SHARED ${SOURCES} ${HEADERS})
set_target_properties(manifoldsTest PROPERTIES
	                  VERSION ${PROJECT_VERSION}
                      SOVERSION 0.1)
#target_compile_definitions(manifoldsTest PRIVATE -DMANIFOLDS_EXPORT -DMNF_ASSERT_THROW -DEIGEN_RUNTIME_NO_MALLOC)
set_target_properties(manifoldsTest PROPERTIES
                    COMPILE_DEFINITIONS "PRIVATE -DMANIFOLDS_EXPORT -DMNF_ASSERT_THROW")

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

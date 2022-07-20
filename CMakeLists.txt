cmake_minimum_required(VERSION 3.5)

project(2dslam LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} /usr/lib/python3.10/site-packages/numpy/core/include)


add_executable(2dslam main.cpp vecext.hpp simulate.hpp preproc.hpp segmentation.hpp diff.hpp dop.hpp ekf.hpp ls.hpp)

# external libs
# find_package(FastAD CONFIG REQUIRED)
find_package(Ceres REQUIRED)
find_package(autodiff REQUIRED)
# matplot
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(Python3 COMPONENTS NumPy)
find_package(OpenMP REQUIRED)

if(Python3_NumPy_FOUND AND OpenMP_CXX_FOUND)
  target_link_libraries(2dslam INTERFACE
    Python3::NumPy OpenMP::OpenMP_CXX FastAD::FastAD autodiff::autodiff
  )
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ")
else()
  target_compile_definitions(2dslam INTERFACE WITHOUT_NUMPY)
endif()
install(
  TARGETS 2dslam
  EXPORT install_targets
)

target_include_directories(2dslam
  PRIVATE
    ${CERES_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(2dslam
  PRIVATE
    ${CERES_LIBRARIES}
    Python3::Python
    Python3::Module
)


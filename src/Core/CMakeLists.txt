add_library(MyFEM_core INTERFACE)
add_library(MyFEM::Core ALIAS MyFEM_core)
target_include_directories(MyFEM_core INTERFACE
${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(MyFEM_core INTERFACE cxx_std_11)
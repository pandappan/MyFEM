add_library(MyFEM_io STATIC
Reader.cpp
Writer.cpp
        VtuExporter.cpp)
add_library(MyFEM::IO ALIAS MyFEM_io)
target_include_directories(MyFEM_io PUBLIC
${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(MyFEM_io PUBLIC
MyFEM::Core
MyFEM::Model)
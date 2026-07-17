add_library(MyFEM_io STATIC
        Writer.cpp
        VtuExporter.cpp
        VtuWriter.cpp
        JsonReader.cpp
)
add_library(MyFEM::IO ALIAS MyFEM_io)
target_include_directories(MyFEM_io PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}
${CMAKE_SOURCE_DIR}/third_party/nlohmann)
target_link_libraries(MyFEM_io PUBLIC
MyFEM::Core
MyFEM::Model)
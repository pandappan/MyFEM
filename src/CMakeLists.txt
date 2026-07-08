add_subdirectory(Analysis)
add_subdirectory(Core)
add_subdirectory(IO)
add_subdirectory(Model)
add_executable(MyFEM main.cpp)
target_link_libraries(MyFEM PRIVATE
        MyFEM::Model
        MyFEM::Analysis
        MyFEM::IO)
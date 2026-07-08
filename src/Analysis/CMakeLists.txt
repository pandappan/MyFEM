add_library(MyFEM_Analysis STATIC
        Assembly.cpp
        Solver.cpp)
add_library(MyFEM::Analysis ALIAS MyFEM_Analysis)
target_include_directories(MyFEM_Analysis PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(MyFEM_Analysis PUBLIC
        MyFEM::Core
        MyFEM::Model)

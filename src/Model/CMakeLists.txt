add_library(MyFEM_model STATIC
        Model.cpp
        Node.cpp
        ElementFactory.cpp

        Element/Element.cpp
        Element/ElementGroup.cpp
        Element/CContinuumElement.cpp
        Element/Bar3D.cpp
        Element/Q4.cpp

        Material/Material.cpp
        Material/BarMaterial.cpp
        Material/CPlaneStressMaterial.cpp
        Material/PlaneStrainMaterial.cpp
        Material/Solid3DMaterial.cpp
)
add_library(MyFEM::Model ALIAS MyFEM_model)

target_include_directories(MyFEM_model PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(MyFEM_model PUBLIC
        MyFEM::Core
)
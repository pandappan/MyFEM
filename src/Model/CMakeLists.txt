add_library(MyFEM_model STATIC
        Model.cpp
        Node.cpp
        ElementFactory.cpp

        Element/Element.cpp
        Element/ElementGroup.cpp
        Element/ContinuumElement.cpp
        Element/Bar3D.cpp
        Element/Q4.cpp
        Element/H8.cpp
        Element/Q4_AX.cpp

        Material/Material.cpp
        Material/BarMaterial.cpp
        Material/PlaneStressMaterial.cpp
        Material/PlaneStrainMaterial.cpp
        Material/Solid3DMaterial.cpp
        Material/AxisymMaterial.cpp
)
add_library(MyFEM::Model ALIAS MyFEM_model)

target_include_directories(MyFEM_model PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/third_party/nlohmann
)
target_link_libraries(MyFEM_model PUBLIC
        MyFEM::Core
)
add_library(Material STATIC
        Material.cpp
        CPlaneStressMaterial.cpp
        BarMaterial.cpp
        PlaneStrainMaterial.cpp
        Solid3DMaterial.cpp)
target_include_directories(Material PUBLIC .)
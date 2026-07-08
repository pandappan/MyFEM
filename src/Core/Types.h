//
// Created by Administrator on 2026/7/5.
//

#pragma once

enum DOFIndex: int {
    UX = 0,
    UY = 1,
    UZ = 2,
    NDOF_MAX=3
};

enum class ElementTypes: int {
    UNDEFINED = 0,
    Bar2D = 1,
    Bar3D = 2,
    Q4_PS = 3,
    Q4_PE = 4,
    T3_PS = 5,
    T3_PE = 6,
    H8 = 7,
    Tet4 =8
};

inline const char* ElementTypeName(ElementTypes type) {
    switch (type) {
        case ElementTypes::Bar2D :
            return "Bar2D";
        case ElementTypes::Bar3D :
            return "Bar3D";
        case ElementTypes::Q4_PS :
            return "Q4_PS";
        case ElementTypes::Q4_PE :
            return "Q4_PE";
        case ElementTypes::T3_PS :
            return "T3_PS";
        case ElementTypes::T3_PE :
            return "T3_PE";
        case ElementTypes::H8 :
            return "H8";
        case ElementTypes::Tet4 :
            return "Tet4";
        default :
            return "undefined element type";
    }
}
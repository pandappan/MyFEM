/**
 * 索引约定
 * 输入文件：单元，节点，面，自由度，全部采用1基，与ABAQUS一致
 * 模型存储：全部采用1基
 * 函数参数：Assembler调用内部函数时，转换为0基
 * eqn[]:方程号0表示非自由，>0表示自由状态
 * DofIndex: 0基
 * 单元节点编号顺序
 * Q4单元
    N4 ────── N3
    │          │
    │          │
    │          │
    N1 ────── N2
    局部节点顺序：N1 → N2 → N3 → N4（逆时针），与ABAQUS一致
        N8──────N7
       /│      /│
      / │     / │
     N5──────N6 │
     │  │    │  │  
     │  N4───│─N3
     │ /     │ /
     N1──────N2
     底面 z=-1：N1(-1,-1,-1) N2(+1,-1,-1) N3(+1,+1,-1) N4(-1,+1,-1)  逆时针
     顶面 z=+1：N5(-1,-1,+1) N6(+1,-1,+1) N7(+1,+1,+1) N8(-1,+1,+1)  逆时针
 * 面-节点顺序
 * Q4单元
                face 2 (N3→N4)
            N4 ─────────────── N3
            │                   │
            │                   │
    face 3  │                   │  face 1
    (N4→N1) │                   │  (N2→N3)，与ABAQUS保持一致
            │                   │
            N1 ─────────────── N2
                  face 0 (N1→N2)
    face 0 : N1-N2-N3-N4 (底 z=-1)
    face 1 : N5-N6-N7-N8 (顶 z=+1)
    face 2 : N1-N2-N6-N5 (前 y=-1)
    face 3 : N2-N3-N7-N6 (右 x=+1)
    face 4 : N3-N4-N8-N7 (后 y=+1)
    face 5 : N4-N1-N5-N8 (左 x=-1)
 * 边界约束约定
 * bcode[0] = 0 : 自由自由度，进入方程
 * bcode[1] = 1 : 固定自由度，不加入方程
 * bcode[2] = 2 : 给定位移自由度，不加入方程，设置给定位移
 **/

#pragma once
#include <string>

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

enum class MaterialTypes: int {
    UNDEFINED = 0,
    Bar = 1,
    PS = 2,
    PE = 3,
    SOLID = 4
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

inline ElementTypes StringToElementType(const std::string& type) {
    if (type == std::string("Bar3D")) return ElementTypes::Bar3D;
    if (type == std::string("Q4_PS")) return ElementTypes::Q4_PS;
    if (type == std::string("Q4_PE")) return ElementTypes::Q4_PE;
    if (type == std::string("H8")) return ElementTypes::H8;
    return ElementTypes::UNDEFINED;
}

inline MaterialTypes StringToMaterialType(const std::string& type) {
    if (type == std::string("bar")) return MaterialTypes::Bar;
    if (type == std::string("plane_stress")) return MaterialTypes::PS;
    if (type == std::string("plane_strain")) return MaterialTypes::PE;
    if (type == std::string("solid3d")) return MaterialTypes::SOLID;
    return MaterialTypes::UNDEFINED;
}

inline bool MaterialCompatibleWithElement(MaterialTypes matType, ElementTypes elemType) {
    switch (elemType) {
        case ElementTypes::Bar3D: return matType == MaterialTypes::Bar;
        case ElementTypes::Q4_PS:
        case ElementTypes::T3_PS: return matType == MaterialTypes::PS;
        case ElementTypes::Q4_PE:
        case ElementTypes::T3_PE: return matType == MaterialTypes::PE;
        case ElementTypes::H8:
        case ElementTypes::Tet4:  return matType == MaterialTypes::SOLID;
        default: return false;
    }
}
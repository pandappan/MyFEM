//
// Created by Administrator on 2026/7/5.
//

#include "ElementFactory.h"
#include "Element/Bar3D.h"
#include "Element/Q4.h"
#include "Element/H8.h"
#include "Material/Material.h"
#include "Material/BarMaterial.h"
#include "Material/CPlaneStressMaterial.h"
#include "Material/PlaneStrainMaterial.h"
#include "Material/Solid3DMaterial.h"
template <class T>
static std::unique_ptr<T> make_unique_() {
    return std::unique_ptr<T>(new T());
}

std::unique_ptr<CElement> CreateElement(ElementTypes type) {
    switch (type) {
        case ElementTypes::Bar3D: return make_unique_<CBar3D>();
        case ElementTypes::Q4_PS: return make_unique_<CQ4>();
        case ElementTypes::Q4_PE: return make_unique_<CQ4>();
        case ElementTypes::H8:    return make_unique_<CH8>();
        default: return nullptr;
    }
}

std::unique_ptr<CMaterial> CreateMaterial(ElementTypes type) {
    switch (type) {
        case ElementTypes::Bar3D: return make_unique_<CBarMaterial>();
        case ElementTypes::Q4_PS: return make_unique_<CPlaneStressMaterial>();
        case ElementTypes::Q4_PE: return make_unique_<CPlaneStrainMaterial>();
        case ElementTypes::H8   : return make_unique_<CSolid3DMaterial>();
        default : return nullptr;
    }
}

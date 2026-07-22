//
// Created by Administrator on 2026/7/5.
//

#include "ElementFactory.h"
#include "Element/Bar3D.h"
#include "Element/Q4.h"
#include "Element/H8.h"
#include "Element/Q4_AX.h"
#include "Material/Material.h"
#include "Material/BarMaterial.h"
#include "Material/CPlaneStressMaterial.h"
#include "Material/PlaneStrainMaterial.h"
#include "Material/Solid3DMaterial.h"
#include "Material/CAxisymMaterial.h"
template <class T>
static std::unique_ptr<T> make_unique_() {
    return std::unique_ptr<T>(new T());
}

std::unique_ptr<CElement> CreateElementByString(const std::string& type) {
    if (type == "Bar3D") {
        return make_unique_<CBar3D>();
    }
    if (type == "Q4_PS" || type == "Q4_PE") {
        return make_unique_<CQ4>();
    }
    if (type == "H8") {
        return make_unique_<CH8>();
    }
    if (type == "Q4_AX") {
        return make_unique_<CQ4_AX>();
    }
    return nullptr;
}

std::unique_ptr<CMaterial> CreateMaterialByString(const std::string &type) {
    if (type == "bar") return make_unique_<CBarMaterial>();
    if (type == "plane_stress") return make_unique_<CPlaneStressMaterial>();
    if (type == "plane_strain") return make_unique_<CPlaneStrainMaterial>();
    if (type == "solid3d") return make_unique_<CSolid3DMaterial>();
    if (type == "axisym") return make_unique_<CAxisymMaterial>();
    return nullptr;
}

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
#include "Material/PlaneStressMaterial.h"
#include "Material/PlaneStrainMaterial.h"
#include "Material/Solid3DMaterial.h"
#include "Material/AxisymMaterial.h"

std::unique_ptr<Element> CreateElementByString(const std::string& type) {
    if (type == "Bar3D") {
        return std::make_unique<Bar3D>();
    }
    if (type == "Q4_PS" || type == "Q4_PE") {
        return std::make_unique<Q4>();
    }
    if (type == "H8") {
        return std::make_unique<H8>();
    }
    if (type == "Q4_AX") {
        return std::make_unique<Q4_AX>();
    }
    return nullptr;
}

std::unique_ptr<Material> CreateMaterialByString(const std::string &type) {
    if (type == "bar") return std::make_unique<BarMaterial>();
    if (type == "plane_stress") return std::make_unique<PlaneStressMaterial>();
    if (type == "plane_strain") return std::make_unique<PlaneStrainMaterial>();
    if (type == "solid3d") return std::make_unique<Solid3DMaterial>();
    if (type == "axisym") return std::make_unique<AxisymMaterial>();
    return nullptr;
}

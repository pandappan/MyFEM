//
// Created by Administrator on 2026/7/15.
//
#include <iostream>
#include <ostream>
#include "JsonReader.h"
#include "Material/Material.h"
#include "Material/BarMaterial.h"
#include "Material/CPlaneStressMaterial.h"
#include "Material/PlaneStrainMaterial.h"
#include "ElementFactory.h"
#include "Element/Element.h"
#include "Model.h"

namespace {
    void FillMaterialParms(CMaterial* mat, const json& matJson) {
        // ! 通用数据
        mat->rho = matJson.value("rho",0.0);
        mat->E = matJson.at("E").get<double>();
        // 杆独有
        if (mat->matType == MaterialTypes::Bar) {
            mat->nu = 0.0;
            dynamic_cast<CBarMaterial*>(mat)->Area = matJson.value("area",1.0);
        } else {
            // 连续介质单元
            mat->nu = matJson.at("nu").get<double>();
            if (mat->matType == MaterialTypes::PS) {
                dynamic_cast<CPlaneStressMaterial*>(mat)->thk = matJson.value("thk",1.0);
            } else if (mat->matType == MaterialTypes::PE) {
                dynamic_cast<CPlaneStrainMaterial*>(mat)->thk = matJson.value("thk",1.0);
            }
        }
    }
}

std::vector<unsigned int> JsonReader::DofStringToInt(const std::vector<std::string>& dofsString) {
    std::vector<unsigned int> dofs_0;
    dofs_0.reserve(dofsString.size());
    for (const auto& dofString : dofsString) {
        if (dofString == "x" || dofString == "X") {
            dofs_0.push_back(0);
        } else if (dofString == "y" || dofString == "Y") {
            dofs_0.push_back(1);
        } else if (dofString == "z" || dofString == "Z") {
            dofs_0.push_back(2);
        } else {
            throw std::runtime_error("DofsStringToInt: Invalid dof string");
        }
    }
    return dofs_0;
}

unsigned int JsonReader::DofStringToInt(const std::string& dofsString) {
    if (dofsString == "x" || dofsString == "X") {
        return 0;
    }
    if (dofsString == "y" || dofsString == "Y") {
        return 1;
    }
    if (dofsString == "z" || dofsString == "Z") {
        return 2;
    }
    throw std::runtime_error("DofStringToInt: Invalid dof string");
}

unsigned int JsonReader::Base1ToBase0(unsigned int idx_1) {
    return idx_1 - 1;
}

std::vector<unsigned int> JsonReader::Base1ToBase0(const std::vector<unsigned int> &idxs_1) {
    std::vector<unsigned int> idxs_0;
    idxs_0.reserve(idxs_1.size());
    for (unsigned int i : idxs_1) {
        idxs_0.push_back(i-1);
    }
    return idxs_0;
}


bool JsonReader::Read(const std::string &filename, Model &model) {
    // 读文件，创建json
    std::ifstream input(filename.c_str());
    if (!input) {
        std::cerr << "Error in Read: " << filename << std::endl;
        return false;
    }
    json j;
    try {
        // 转为json
         j = json::parse(input);
        // 将json读取文件中
        if (!ParseHeader(j,model)) {return false;}
        if (!ParseNodes(j,model)) {return false;}
        if (!ParseMaterials(j,model)) {return false;}
        if (!ParseElementGroups(j,model)) {return false;}
        if (!ParseBoundaryConditions(j,model)) {return false;}
        if (!ParseLoads(j,model)) {return false;}
    } catch (const json::parse_error& e) {
        std::cerr << "Json parse error: " << e.what() << std::endl;
        return false;
    } catch (const std::exception& e) {
        std::cerr << "Input error: " << e.what() << std::endl;
        return false;
    }
    return true;
}

bool JsonReader::ParseHeader(const json &j, Model &model) {
    model.title = j.at("title").get<std::string>();
    model.dimension = j.at("dimension").get<unsigned int>();
    if (model.dimension !=2 && model.dimension !=3) {
        throw std::runtime_error("Error: Only 2 or 3 dimensions are supported");
    }
    model.modex = 1;
    return true;
}

// 设置节点的约束情况
// 约定：节点编号连续
bool JsonReader::ParseNodes(const json &j, Model &model) {
    const auto& nodesJson = j.at("nodes");
    model.nodes.clear();
    model.nodes.reserve(nodesJson.size());
    std::vector<double> XYZ(3);
    unsigned int globalNodeId = 0;
    for (auto& nodeJson : nodesJson) {
        unsigned int nodeId_0 = Base1ToBase0(nodeJson.at("id").get<unsigned int>());
        // 检测节点是否连续
        if (nodeId_0 != globalNodeId++) {
            throw std::runtime_error("ParseNodes: node id must be contiguous" + std::to_string(globalNodeId));
        }
        XYZ[0] = nodeJson.at("x").get<double>();
        XYZ[1] = nodeJson.at("y").get<double>();
        XYZ[2] = nodeJson.value("z", 0.0);
        CNode node;
        // 设置节点基本信息
        node.SetGeom(nodeId_0, XYZ);
        // 默认全部约束保持自由，根据维度锁定Z自由度
        node.SetDimConstraints(model.dimension);
        model.nodes.push_back(node);
    }
    return true;
}

bool JsonReader::ParseMaterials(const json& j, Model& model) {
    // 初始化容器
    model.materials.clear();
    model.materials.reserve(j.at("materials").size());
    // 读取材料数据
    for (const auto& m: j.at("materials")) {
        std::string matString = m.at("type");
        auto mat = CreateMaterialByString(matString);
        if (!mat) {
            throw std::runtime_error("ParseMaterials: Material \"" + matString + "\" not found");
        }
        FillMaterialParms(mat.get(), m);
        model.materials.push_back(std::move(mat));
    }
    return true;
}


/*
 *elemeng_groups:[
 *{type:
 *material:
 *element:[{id:, connectivity:[]},
 *{id:, connectivity:[]}]}
 *]
 */
// 全局编号-单元指针映射表，初始化单元基础信息
bool JsonReader::ParseElementGroups(const json &j, Model &model) {
    const auto& groupsJson = j.at("element_groups");
    // 1.初始化单元组数
    // 2.统计总单元数，初始化映射表
    unsigned int totalGroups = 0;
    totalGroups = groupsJson.size();
    model.groups.clear();
    model.groups.resize(totalGroups);
    // 初始化组内单元信息，添加映射表数据
    unsigned int idxGroup = 0;
    for (const auto& g : groupsJson) {
        // 设置组内的信息
        std::string elemString = g.at("type").get<std::string>();
        ElementTypes elemType = StringToElementType(elemString);
        if (elemType == ElementTypes::UNDEFINED) {
            throw std::runtime_error("ParseElementGroups: unknow element type: " + elemString);
        }
        // 组内的材料信息
        unsigned int matId_0 = Base1ToBase0(g.at("material").get<unsigned int>());
        if (matId_0 >= model.materials.size()) {
            throw std::runtime_error("ParseElementGroups: Invalid material reference " + std::to_string(matId_0));
        }
        if (!MaterialCompatibleWithElement(model.materials[matId_0]->matType, elemType))
            throw std::runtime_error("Material incompatible with element type");
        CMaterial* matPtr = model.materials[matId_0].get();
        CElementGroup& group = model.groups[idxGroup];
        group.SetGroupsInfo(elemType,g.at("elements").size());
        // 初始化组内单元信息
        for (const auto& e: g.at("elements")) {
            auto elemId_0 = Base1ToBase0(e.at("id").get<unsigned int>());
            auto connectivity_0 =
                Base1ToBase0(e.at("connectivity").get<std::vector<unsigned int>>());
            auto elem = CreateElementByString(elemString);
            if (!elem) {
                throw std::runtime_error("Error: Could not create element \"" + elemString + "\"");
            }
            // 单元基础信息，单元积分点信息
            elem->SetElementType(elemType);
            elem->SetElementInfo(elemId_0, matPtr, connectivity_0,model.nodes);
            elem->OnSetupComplete();
            // 加入单元组
            group.AddElement(std::move(elem));
        }
        idxGroup++;
    }
    return true;
}

bool JsonReader::ParseBoundaryConditions(const json &j, Model &model) {
    // 整体检查
    if (!j.contains("boundary_conditions")) {
        throw std::runtime_error("Error: No boundary conditions provided");
    }
    const auto& BCsJson = j.at("boundary_conditions");
    // 固定约束，节点编号，自由度固定约束情况
    if (BCsJson.contains("fixed")) {
        for (const auto& item: BCsJson.at("fixed")) {
            const auto& nodesId_1 = item.at("nodes").get<std::vector<unsigned int>>();
            const auto& dofsString = item.at("dof").get<std::vector<std::string>>();
            const auto& dofs_0 = DofStringToInt(dofsString);
            for (auto nodeId_1 : nodesId_1) {
                unsigned int nodeId_0 = nodeId_1 - 1;
                model.nodes[nodeId_0].SetFixConstraints(dofs_0);
            }
        }
    }
    // 指定位移约束
    if (BCsJson.contains("prescribed")) {
        for (const auto& item: BCsJson.at("prescribed")) {
            const auto& nodesId_0 = Base1ToBase0(item.at("node").get<unsigned int>());
            const auto& dofString = item.at("dof").get<std::string>();
            const auto& dof_0 = DofStringToInt(dofString);
            double value  = item.at("value").get<double>();
            model.nodes[nodesId_0].SetPreDispConstraints(dof_0, value);
        }
    }
    return true;
}

bool JsonReader::ParseLoads(const json &j, Model &model) {
    // 不包含载荷
    if (!j.contains("loads")) return true;
    const auto& loadsJson = j.at("loads");
    // 点载荷，直接进入节点中
    if (loadsJson.contains("concentrated")) {
        model.cloads.clear();
        model.cloads.reserve(loadsJson.at("concentrated").size());
        for (auto& cloadJson: loadsJson.at("concentrated")) {
            unsigned int nodeId_0 = Base1ToBase0(cloadJson.at("node").get<unsigned int>());
            // 检查节点索引是否超过范围
            if (nodeId_0 >= model.nodes.size()) {
                throw std::runtime_error("ParseLoads: node id is out of range " + std::to_string(nodeId_0+1));
            }
            unsigned int dof_0 = DofStringToInt(cloadJson.at("dof").get<std::string>());
            double value = cloadJson.at("value").get<double>();
            model.nodes[nodeId_0].AddForce(dof_0, value);
        }
    }
    // 面/线载荷
    if (loadsJson.contains("surface")) {
        const auto& surfacesJson = loadsJson.at("surface");
        model.sloads.clear();
        model.sloads.reserve(surfacesJson.size());
        for (auto& surfaceJson : surfacesJson) {
            SurfaceLoad sload{};
            sload.elemId_0 = Base1ToBase0(surfaceJson.at("element").get<int>());
            // 检测映射表中是否包含该单元编号
            auto it = model.globalElementMap.find(sload.elemId_0);
            if (it == model.globalElementMap.end()) {
                throw std::runtime_error("ParseLoads: surface does not find elemId: " + std::to_string(sload.elemId_0));
            }
            sload.faceId_0 = Base1ToBase0(surfaceJson.at("face").get<int>());
            sload.dof_0 = DofStringToInt(surfaceJson.at("dof").get<std::string>());
            sload.value = surfaceJson.at("value").get<double>();
            model.sloads.push_back(sload);
        }
    }
    // 体载荷
    if (loadsJson.contains("body_force")) {
        const auto& bodyForceJson = loadsJson.at("body_force");
        for (unsigned int i = 0; i < bodyForceJson.size(); i++) {
            model.bodyForce[i] = bodyForceJson.at(i).get<double>();
        }
    }
    return true;
}
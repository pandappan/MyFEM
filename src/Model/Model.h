//
// Created by Administrator on 2026/7/5.
//

#pragma once
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include "Node.h"
#include "Material/Material.h"
#include "Element/ElementGroup.h"
#include "../Core/SkylineMatrix.h"


struct ConcentratedLoad {
    unsigned int node_0;
    unsigned int dof_0;
    double value;
};

struct PreDisplacement {
    unsigned int node_0;
    unsigned int dof_0;
    double value;
};

struct SurfaceLoad {
    unsigned int elemId_0;
    unsigned int faceId_0;
    unsigned int dof_0;
    double       value;   ///< Distributed load per unit length/area
};

struct MPCTerm {
    unsigned int node_0;
    unsigned int dof_0;
    double coeff;
};

struct MPC {
    unsigned int slaveNode_0;
    unsigned int slaveDof_0;
    std::vector<MPCTerm> masters;
    double beta = 0.0;
};

class Model {
public:
    // 模型网格，载荷数据
    std::string title;
    unsigned int dimension = 3;
    unsigned int modex = 0;
    std::vector<std::unique_ptr<Material>> materials;
    std::vector<Node> nodes;
    std::vector<CElementGroup> groups;
    std::vector<ConcentratedLoad> cloads;
    std::vector<PreDisplacement> predisplacements;
    std::vector<SurfaceLoad> sloads;
    std::unordered_map<unsigned int, Element*> globalElementMap; // 单元全局编号(0基)-单元指针映射
    double bodyForce[3] = {0.0, 0.0, 0.0};
    // 多点约束方程
    std::vector<MPC> mpcs;
    // 从节点的自由度-全局映射
    std::unordered_map<unsigned int, unsigned int> slaveDofToMpc;
    // 整体刚度矩阵，右端项
    unsigned int neq = 0;
    std::unique_ptr<SkylineMatrix<double>> K;
    std::vector<double> force; // 求解前为外载荷，求解后存储节点位移
    // 查询函数
    unsigned int GetNumNodes() const {return nodes.size();}
    unsigned int GetNumGroups() const {return groups.size();}
    Material* GetMaterialPtr(unsigned int index0) const {return materials[index0].get();}
    int FindMPCBySlave(unsigned int slaveNode_0, unsigned int slaveDof_0) const;
    Model() = default;
    ~Model();
    // 不允许拷贝，只允许移动
    Model(const Model&) = delete;
    Model& operator=(const Model&) = delete;
    Model(Model&&) = default;
    Model& operator=(Model&&) = default;
};

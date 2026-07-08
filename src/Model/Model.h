//
// Created by Administrator on 2026/7/5.
//

#pragma once
#include <memory>
#include <string>
#include <vector>
#include "Node.h"
#include "Element/ElementGroup.h"
#include "Material/Material.h"
#include "../Core/SkylineMatrix.h"

struct ConcentratedLoad {
    unsigned int node;
    unsigned int dof;
    double value;
};

struct PreDisplacement {
    unsigned int node;
    unsigned int dof;
    double value;
};

class Model {
public:
    // 模型网格，载荷数据
    std::string title;
    unsigned int dimension = 3;
    unsigned int modex = 0;
    std::vector<CNode> nodes;
    std::vector<CElementGroup> groups;
    std::vector<ConcentratedLoad> cloads;
    std::vector<PreDisplacement> predisplacements;
    // 整体刚度矩阵，右端项
    unsigned int neq = 0;
    std::unique_ptr<CSkylineMatrix<double>> K;
    std::vector<double> force; // 求解前为外载荷，求解后存储节点位移
    // 查询函数
    unsigned int GetNumNodes() const {return nodes.size();}
    unsigned int GetNumGroups() const {return groups.size();}
    Model() = default;
    ~Model();
    // 不允许拷贝，只允许移动
    Model(const Model&) = delete;
    Model& operator=(const Model&) = delete;
    Model(Model&&) = default;
    Model& operator=(Model&&) = default;
};

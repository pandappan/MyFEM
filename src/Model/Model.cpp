//
// Created by Administrator on 2026/7/5.
//
#include "Model.h"
#include "Material/Material.h"
#include "../Core/SkylineMatrix.h"


int Model::FindMPCBySlave(unsigned int slaveNode_0, unsigned int slaveDof_0) {
    // 查询MPC的表，如果能查到则正常输出
    unsigned int key = slaveDof_0 * CNode::NDF + slaveDof_0;
    auto it = slaveDofToMpc.find(key);
    if (it == slaveDofToMpc.end()) {
        return -1;
    } else {
        return static_cast<int>(it->second);
    }
}

Model::~Model() = default;

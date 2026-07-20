//
// Created by Administrator on 2026/6/11.
//
#include <iostream>
#include "Element.h"
#include "../Node.h"
#include "../../Core/Types.h"
#include "../../Core/DenseMatrix.h"

// 初始化单元基础信息
void CElement::SetElementInfo(unsigned int elemId_0, CMaterial* matPtr,
    const std::vector<unsigned int>& connectivity_0, std::vector<CNode>& nodeList) {
    ElementNumber_   = elemId_0;
    ElementMaterial_ = matPtr;
    for (unsigned int i = 0; i < NEN_; ++i)
        nodes_[i] = &nodeList[connectivity_0[i]];
    // 激活单元所述节点的自由度
    RegisterDofsOnNodes();
}


// 获取单元所有节点的坐标
// x_iI
CElement::CElement():elementType_(ElementTypes::UNDEFINED),ElementNumber_(0),
NDim_(0), NEN_(0), ND_(0), ElementMaterial_(nullptr), volume_((0.0)){}

DenseMatrix<double> CElement::GetNodeCoordinates() const {
    DenseMatrix<double> nodeCoords(NDim_, NEN_);
    for (unsigned int k = 0; k < NEN_; k++) {
        for (unsigned int d = 0; d < NDim_; d++) {
            nodeCoords(d, k) = nodes_[k]->XYZ[d];
        }
    }
    return nodeCoords;
}


// 根据单元-节点连接关系，形成单元的定位数组，用于形成总体刚度矩阵
void CElement::GenerateLocationMatrix()
{
    const DOFIndex* activeDOFs  = GetActiveDOFs();
    unsigned int ndofs = GetNumActiveDOFsPerNode();
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < ndofs; D++)
            LocationMatrix_[i++] = nodes_[N]->eqn[activeDOFs[D]];
}

void CElement::AllocateStorage(unsigned int nDim, unsigned int nen, unsigned int nd) {
    NEN_ = nen;
    NDim_ = nDim;
    ND_ = nd;
    nodes_.assign(NEN_, nullptr);
    LocationMatrix_.assign(ND_, 0);
}

void CElement::GetElementNodesDisp(std::vector<double>& nodesDisp,
                                   std::vector<int>& nodesBcode) const {
    const DOFIndex* dofs = GetActiveDOFs();
    unsigned int    ndof = GetNumActiveDOFsPerNode();
    unsigned int index = 0;
    for (const auto* node : nodes_) {
        for (unsigned int d = 0; d < ndof; ++d) {
            nodesDisp[index]  = node->displacement[dofs[d]];
            nodesBcode[index] = static_cast<int>(node->bcode[dofs[d]]);
            index++;
        }
    }
}

void CElement::GetElementNodesForce(std::vector<double>& nodesForce) {
    const DOFIndex* dofs = GetActiveDOFs();
    unsigned int    ndof = GetNumActiveDOFsPerNode();
    unsigned int index = 0;
    for (const auto* node : nodes_) {
        for (unsigned int d = 0; d < ndof; ++d) {
            nodesForce[index] = node->nodeForce[dofs[d]];
            index++;
        }
    }
}

// 反力：R_i = Σ_j K(i,j)·u_j - f_i
// 只写回固定、指定位移约束的反力，自由，从属位移则不处理
void CElement::CalculateBCForce() {
    DenseMatrix<double> Ke(ND_, ND_);
    Ke.SetZero();
    ElementStiffness(Ke);

    std::vector<double> nodesDisp(ND_);
    std::vector<int>    nodesBcode(ND_);
    GetElementNodesDisp(nodesDisp, nodesBcode);

    std::vector<double> nodesForce(ND_);
    GetElementNodesForce(nodesForce);

    const DOFIndex*    dofs = GetActiveDOFs();
    const unsigned int ndof = GetNumActiveDOFsPerNode();

    for (unsigned int i = 0; i < ND_; ++i) {
        // 只对固定和指定位移求解反力，自由和从属位移则跳过
        if (nodesBcode[i] != 1 && nodesBcode[i] != 2) continue;

        double KU = 0.0;
        for (unsigned int j = 0; j < ND_; ++j)
            KU += Ke(i, j) * nodesDisp[j];

        double BCForce = KU - nodesForce[i];

        unsigned int nodeIdx  = i / ndof;
        unsigned int localDof = i % ndof;
        nodes_[nodeIdx]->AddBcForce(dofs[localDof], BCForce);
    }
}

bool CElement::CalculateSurfaceLoad(unsigned int faceID, unsigned int dof, double value) {
    if (elementType_ == ElementTypes::Bar3D) {
        std::cerr << "Bar3D do not have surface load" << std::endl;
        return false;
    }
    return true;
}

// 空实现
void CElement::CalculateBodyForce(const double* bodyForce) {
    return;
}

void CElement::SetElementType(ElementTypes elementType) {
    elementType_ = elementType;
}

// 计算单元内的应变能 0.5 u^t K u
double CElement::CalculateElementEnergy() const {
    // 单元所有节点的位移
    std::vector<double> nodesDisp(ND_);
    std::vector<int> nodesBcode(ND_);
    GetElementNodesDisp(nodesDisp, nodesBcode);
    // 单元刚度矩阵
    DenseMatrix<double> Ke(ND_, ND_);
    ElementStiffness(Ke);
    // 调用应变能公式
    double energy = 0.0;
    for (unsigned int i = 0; i < ND_; ++i) {
        for (unsigned int j = 0; j < ND_; ++j) {
            energy += 0.5 * nodesDisp[i] * Ke(i,j) * nodesDisp[j];
        }
    }
    return energy;
}

void CElement::SetupForTesting(std::vector<CNode*> NodeList, CMaterial* Material_) {
    nodes_ = std::move(NodeList);
    ElementMaterial_ = Material_;
}

void CElement::RegisterDofsOnNodes() {
    const DOFIndex* dofs = GetActiveDOFs();
    const unsigned int ndof = GetNumActiveDOFsPerNode();
    for (CNode* node: nodes_) {
        for (unsigned int d = 0; d < ndof; ++d) {
            node->ActivateDof(static_cast<unsigned int>(dofs[d]));
        }
    }
}

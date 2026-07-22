//
// Created by Administrator on 2026/6/10.
//

#include <iostream>
#include <iomanip>
#include "Q4.h"
#include "../Node.h"
#include "../Material/Material.h"
#include "CContinuumElement.h"
#include "ElementGroup.h"

const DOFIndex CQ4::ActiveDOFs[2] = {UX, UY};
const unsigned int CQ4::NumActiveDOFsPerNode = 2;
// 构造函数
CQ4::CQ4() {
    AllocateStorage(2,4,8);
}

// Q4作为面元时，只需要节点数据便可进行初始化
// 只传入节点指针数组，材料数组保持为空，无需生成方程号，初始化积分点处数据，供等效面力使用
void CQ4::AsFaceElem(const std::vector<CNode*> &nodelist) {
    for (int i = 0; i < 4; i++) {
        nodes_[i] = nodelist[i];
    }
    InitializeIntegrationPoints();
}

void CQ4::Write(std::ostream& Output) const {
    Output << std::setw(6) << ElementNumber_;
    for (unsigned int i = 0; i < 4; i++) {
        Output << std::setw(6) << nodes_[i]->Index;
    }
    Output << std::setw(6) << ElementMaterial_->nset << std::endl;
}

void CQ4::WriteElementStress(std::ostream& out) const {
    for (unsigned int ip = 0; ip < GetNumIntegrationPoints(); ++ip) {
        auto stress = ComputeStressAtIntegrationPoint(ip);
        out << std::setw(6) << ElementNumber_
            << std::setw(6) << ip + 1;
        for (double s : stress) out << std::setw(6) << s;
        out << std::endl;
    }
}

GaussData CQ4::GetIntegrationPoint() const {
    GaussData g;
    g.GaussWeights.assign(4,1.0);
    g.GaussPoints.Resize(4,2);
    const double a = 0.577350269189626;
    g.GaussPoints(0, 0) = -a;  g.GaussPoints(0, 1) = -a;
    g.GaussPoints(1, 0) =  a;  g.GaussPoints(1, 1) = -a;
    g.GaussPoints(2, 0) =  a;  g.GaussPoints(2, 1) =  a;
    g.GaussPoints(3, 0) = -a;  g.GaussPoints(3, 1) =  a;
    return g;
}

void CQ4::ComputeShapeFunctions(const std::vector<double>& xi,
    std::vector<double>& N) const {
    N[0] = 0.25 * (1.0 - xi[0]) * (1.0 - xi[1]);
    N[1] = 0.25 * (1.0 + xi[0]) * (1.0 - xi[1]);
    N[2] = 0.25 * (1.0 + xi[0]) * (1.0 + xi[1]);
    N[3] = 0.25 * (1.0 - xi[0]) * (1.0 + xi[1]);
}

void CQ4::ComputeShapeDerivatives(const std::vector<double>& xi,
    DenseMatrix<double>& dN_dxi) const {
    // dN_I/dxi_i
    dN_dxi(0,0) = - 0.25 * (1.0 - xi[1]);
    dN_dxi(0,1) =   0.25 * (1.0 - xi[1]);
    dN_dxi(0,2) =   0.25 * (1.0 + xi[1]);
    dN_dxi(0,3) = - 0.25 * (1.0 + xi[1]);

    dN_dxi(1,0) = - 0.25 * (1.0 - xi[0]);
    dN_dxi(1,1) = - 0.25 * (1.0 + xi[0]);
    dN_dxi(1,2) =   0.25 * (1.0 + xi[0]);
    dN_dxi(1,3) =   0.25 * (1.0 - xi[0]);
}

// faceID: 0基
bool CQ4::CalculateSurfaceLoad(unsigned int faceID, unsigned int dof, double value) {
    // 面局部节点
    std::vector<int> nodesLocalID = GetFaceNodesLocalID(faceID);
    double x1 = nodes_[nodesLocalID[0]]->XYZ[0];
    double y1 = nodes_[nodesLocalID[0]]->XYZ[1];
    double x2 = nodes_[nodesLocalID[1]]->XYZ[0];
    double y2 = nodes_[nodesLocalID[1]]->XYZ[1];
    double len = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    // 积分方案
    std::vector<double> xi(2);
    xi[0] = -0.57735027;
    xi[1] = 0.57735027;
    double weight = 1.0;
    DenseMatrix<double> N(2,2); // (ng, np)
    N(0,0) = (1.0 - xi[0]) * 0.5;
    N(0,1) = (1.0 + xi[0]) * 0.5;
    N(1,0) = (1.0 - xi[1]) * 0.5;
    N(1,1) = (1.0 + xi[1]) * 0.5;
    // 面力等效节点力
    for (unsigned int i = 0; i < 2; i++) {
        double eqforce = 0.5 * len * value * (N(0, i) +N(1,i)) * weight;
        // 等效节点力写入节点中
        nodes_[nodesLocalID[i]]->AddForce(dof, eqforce);
    }
    return true;
}

// FaceID: 0基
std::vector<int> CQ4::GetFaceNodesLocalID(unsigned int faceID) const {
    // Q4 face convention (0-based, CCW):
    //   face 0: N1(0) -> N2(1)  (bottom)
    //   face 1: N2(1) -> N3(2)  (right)
    //   face 2: N3(2) -> N4(3)  (top)
    //   face 3: N4(3) -> N1(0)  (left, closing edge)
    std::vector<int> nodesLocalID(2);
    switch (faceID) {
        case 0: nodesLocalID={0,1}; break;
        case 1: nodesLocalID={1,2}; break;
        case 2: nodesLocalID={2,3}; break;
        case 3: nodesLocalID={3,0}; break;
        default:
            throw std::out_of_range("CQ4::GetFaceNodesLocalID invalid faceID");
    }
    return nodesLocalID;
}

void CQ4::GetVisualizationNodes(DenseMatrix<double>& coords) const {
    coords.Resize(2,NEN_);
    for (unsigned int i = 0; i < NEN_; i++) {
        coords(0,i) = nodes_[i]->XYZ[0];
        coords(1,i) = nodes_[i]->XYZ[1];
    }
}

//! 返回节点变形坐标
void CQ4::GetVisualizationDeformeNodes(DenseMatrix<double>& deformeCoords) const {
    deformeCoords.Resize(2, NEN_);
    for (unsigned int i = 0; i < NEN_; i++) {
        deformeCoords(0,i) = nodes_[i]->XYZ[0] + nodes_[i]->displacement[0];
        deformeCoords(1,i) = nodes_[i]->XYZ[1] + nodes_[i]->displacement[1];
    }
}

//! 单元外推矩阵
DenseMatrix<double> CQ4::GetExprapolationMatrix() const {
    DenseMatrix<double> E(4,4);
    const double sqrt3 = 1.7320508075688772;
    const double a = 1.0 + sqrt3 / 2.0;      // ≈ 1.866
    const double b = -0.5;                   // = -0.5
    const double c = 1.0 - sqrt3 / 2.0;      // ≈ 0.134
    E(0, 0) = a;  E(0, 1) = b;  E(0, 2) = c;  E(0, 3) = b;
    E(1, 0) = b;  E(1, 1) = a;  E(1, 2) = b;  E(1, 3) = c;
    E(2, 0) = c;  E(2, 1) = b;  E(2, 2) = a;  E(2, 3) = b;
    E(3, 0) = b;  E(3, 1) = c;  E(3, 2) = b;  E(3, 3) = a;
    return E;
}

MaterialCategory CQ4::GetRequiredMaterial() const {
    if (elementType_ == ElementTypes::Q4_PE) return MaterialCategory::MechanicalPlaneStrain;
    if (elementType_ == ElementTypes::Q4_PS) return MaterialCategory::MechanicalPlaneStress;
    return MaterialCategory::UNDEFINED;
}

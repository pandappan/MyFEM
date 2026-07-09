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

bool CQ4::Read(std::ifstream& Input, CElementGroup& group, std::vector<CNode>& nodelist) {
    unsigned int MSet;
    std::vector<unsigned int> NodeNum(4);
    for (int i = 0; i < 4; i++) {
        Input >> NodeNum[i];
    }
    Input >> MSet;
    ElementMaterial_ = &(group.GetMaterial(MSet - 1));
    for (int i = 0; i < 4; i++) {
        nodes_[i] = &nodelist[NodeNum[i]-1];
    }
    InitializeIntegrationPoints();
    return true;
}

void CQ4::Write(std::ostream& Output) const {
    Output << std::setw(6) << ElementNumber_;
    for (unsigned int i = 0; i < 4; i++) {
        Output << std::setw(6) << nodes_[i]->NodeNumber;
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

std::vector<int> CQ4::GetFaceNodesLocalID(unsigned int faceID) {
    std::vector<int> nodesLocalID(2);
    switch (faceID) {
        case 0: nodesLocalID={0,1}; break;
        case 1: nodesLocalID={1,2}; break;
        case 2: nodesLocalID={2,3}; break;
        case 3: nodesLocalID={3,1}; break;
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
        deformeCoords(0,i) = nodes_[i]->XYZ[0] + nodes_[i]->Displacement[0];
        deformeCoords(1,i) = nodes_[i]->XYZ[1] + nodes_[i]->Displacement[1];
    }
}
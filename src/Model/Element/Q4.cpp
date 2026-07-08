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
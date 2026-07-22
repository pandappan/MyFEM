//
// Created by Administrator on 2026/7/22.
//

#include "Q4_AX.h"

#include "Node.h"

void CQ4_AX::ComputeBMatrix(unsigned int ip, DenseMatrix<double> &B) const {
    B.SetZero();
    double radius = GetRadiusAtIntegrationPoint(ip);
    const auto& N = integrationPoints_[ip].N;
    const auto& dN_dx = integrationPoints_[ip].dN_dx;
    for (unsigned int I = 0; I < NEN_; I++) {
        double n = N[I];
        double dNdx = dN_dx(0,I);
        double dNdy = dN_dx(1, I);
        B(0, 2*I    ) = dNdx;
        B(1, 2*I + 1) = dNdy;
        B(2, 2*I    ) = n / radius;
        B(3, 2*I    ) = dNdy;
        B(3, 2*I + 1) = dNdx;
    }
}

double CQ4_AX::GetIntegrationVolumeFactor(unsigned int ip) const {
    return 2.0 * 3.14159265358979323846 * GetRadiusAtIntegrationPoint(ip);
}

// faceID: 0基
// 面载荷：约定 value 是单位面积 traction（不含 2π·r）
bool CQ4_AX::CalculateSurfaceLoad(unsigned int faceID, unsigned int dof, double value) {
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
    std::vector<double> weights = {
        2.0 * 3.14159265358979323846 * GetRadiusAtIntegrationPoint(0),
        2.0 * 3.14159265358979323846 * GetRadiusAtIntegrationPoint(1)
    };
    DenseMatrix<double> N(2,2); // (ng, np)
    N(0,0) = (1.0 - xi[0]) * 0.5;
    N(0,1) = (1.0 + xi[0]) * 0.5;
    N(1,0) = (1.0 - xi[1]) * 0.5;
    N(1,1) = (1.0 + xi[1]) * 0.5;
    // 面力等效节点力
    for (unsigned int i = 0; i < 2; i++) {
        double eqforce = 0.5 * len * value * (N(0, i) +N(1,i)) * weights[i];
        // 等效节点力写入节点中
        nodes_[nodesLocalID[i]]->AddForce(dof, eqforce);
    }
    return true;
}

MaterialCategory CQ4_AX::GetRequiredMaterial() const {
    return MaterialCategory::MechanicalAxisym;
}

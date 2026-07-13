#include <iostream>
#include <ostream>
#include <fstream>
#include "H8.h"
#include "Node.h"
#include "ElementGroup.h"
#include "Material/Material.h"
#include "Q4.h"

const DOFIndex CH8::ActiveDOFs[3] = {UX,UY,UZ};
const unsigned int CH8::NumActiveDOFsPerNode = 3;
const int CH8::nodeSign_[8][3] = {
    {-1, -1, -1},
    { 1, -1, -1},
    { 1,  1, -1},
    {-1,  1, -1},
    {-1, -1,  1},
    { 1, -1,  1},
    { 1,  1,  1},
    {-1,  1,  1},
};

CH8::CH8() {
    AllocateStorage(3, 8, 24);
}

bool CH8::Read(std::ifstream &Input, CElementGroup &group, std::vector<CNode> &nodelist) {
    unsigned int MSet;
    std::vector<unsigned int> NodeNum(8);
    for (int i = 0; i < 8; i++) {
        Input >> NodeNum[i];
    }
    Input >> MSet;
    ElementMaterial_ = &(group.GetMaterial(MSet - 1));
    for (int i = 0; i < 8; i++) {
        nodes_[i] = &nodelist[NodeNum[i]-1];
    }
    InitializeIntegrationPoints();
    return true;
}

void CH8::Write(std::ostream &output) const {
    output << std::setw(6) << ElementNumber_;
    for (unsigned int i = 0; i < 8; i++) {
        output << std::setw(6) << nodes_[i]->NodeNumber;
    }
    output << std::setw(6) << ElementMaterial_->nset << std::endl;
}

void CH8::WriteElementStress(std::ostream& out) const {
    for (unsigned int ip = 0; ip < GetNumIntegrationPoints(); ++ip) {
        auto stress = ComputeStressAtIntegrationPoint(ip);
        out << std::setw(6) << ElementNumber_
            << std::setw(6) << ip + 1;
        for (double s : stress) out << std::setw(6) << s;
        out << std::endl;
    }
}

// 所有的高斯积分点以及权重
GaussData CH8::GetIntegrationPoint() const {
    GaussData g;
    g.GaussWeights.assign(8,1.0);
    g.GaussPoints.Resize(8,3); // (nGp, nDim)
    const double a = 0.577350269189626;
    for (unsigned int i = 0; i < 8; i++) {
        for (unsigned int d = 0; d < 3; d++) {
            g.GaussPoints(i, d) = a * nodeSign_[i][d];
        }
    }
    return g;
}

void CH8::ComputeShapeFunctions(const std::vector<double>& xi,
    std::vector<double>& N) const {
    for (unsigned int i = 0; i < 8; i++) {
        const double sx = nodeSign_[i][0];
        const double sy = nodeSign_[i][1];
        const double sz = nodeSign_[i][2];
        N[i] = 0.125 * (1.0 + sx * xi[0]) * (1.0 + sy * xi[1]) * (1.0 + sz * xi[2]);
    }
}

void CH8::ComputeShapeDerivatives(const std::vector<double>& xi,
    DenseMatrix<double>& dN_dxi) const {
    // dN_I/dxi_i
    for (unsigned int i = 0; i < 8; i++) {
        const double sx = nodeSign_[i][0];
        const double sy = nodeSign_[i][1];
        const double sz = nodeSign_[i][2];
        dN_dxi(0,i) = 0.125 * sx * (1.0 + sy * xi[1]) * (1.0 + sz * xi[2]);
        dN_dxi(1,i) = 0.125 * sy * (1.0 + sx * xi[0]) * (1.0 + sz * xi[2]);
        dN_dxi(2,i) = 0.125 * sz * (1.0 + sx * xi[0]) * (1.0 + sy * xi[1]);
    }
}

// 面载功能先不加入

void CH8::GetVisualizationNodes(DenseMatrix<double>& coords) const {
    coords.Resize(3,NEN_);
    for (unsigned int i = 0; i < NEN_; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            coords(j,i) = nodes_[i]->XYZ[j];
        }
    }
}

void CH8::GetVisualizationDeformeNodes(DenseMatrix<double>& deformeCoords) const {
    deformeCoords.Resize(3, NEN_);
    for (unsigned int i = 0; i < NEN_; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            deformeCoords(j,i) = nodes_[i]->XYZ[j] + nodes_[i]->Displacement[j];
        }
    }
}

// 面载转化为等效节点载荷
bool CH8::CalculateSurfaceLoad(unsigned int faceID, unsigned int dof, double value) {
    // 从局部面元编号转换为节点数组
    std::vector<int> nodesLocalID = GetFaceNodesLocalID(faceID);
    std::vector<CNode*> nodeList(4);
    for (unsigned int i = 0; i < 4; i++) {
        nodeList[i] = nodes_[i];
    }
    // 构建局部面元，并初始化形函数等信息
    CQ4 faceElem;
    faceElem.AsFaceElem(nodeList);
    // 计算等效节点载荷并输出到节点
    const auto& gsData =faceElem.GetIntegrationPoints();
    for (unsigned int n = 0; n < 4; n++) {
        double eqforce = 0.0;
        for (unsigned int ip = 0; ip < 4; ip++) {
                eqforce += gsData[ip].N[n] * value * gsData[ip].detJ_times_weight;
        }
        // 通过当前H8单元输入节点力中
        nodes_[nodesLocalID[n]]->AddForce(dof, eqforce);
    }
    return true;
}

// H8单元面的局部节点编号
std::vector<int> CH8::GetFaceNodesLocalID(unsigned int faceID) const {
    std::vector<int> nodesLocalID(4);
    switch (faceID) {
        case 0: nodesLocalID = {0, 1, 2, 3}; break;
        case 1: nodesLocalID = {5, 6, 7, 8}; break;
        case 2: nodesLocalID = {0, 1, 5, 4}; break;
        case 3: nodesLocalID = {1, 2, 6, 5}; break;
        case 4: nodesLocalID = {2, 3, 7, 6}; break;
        case 5: nodesLocalID = {3, 0, 4, 7}; break;
        default: {
            throw std::out_of_range("GetFaceNodesLocalID: Unknown faceID");
        }
    }
    return nodesLocalID;
}

//! 单元外推矩阵
DenseMatrix<double> CH8::GetExprapolationMatrix() const {
    DenseMatrix<double> E(8,8);
    const double sqrt3 = 1.7320508075688772;
    const double a = 0.5 * (1.0 + sqrt3);
    const double b = 0.5 * (1.0 - sqrt3);
    for (unsigned int i = 0; i < 8; i++) {
        for (unsigned int gp = 0; gp < 8; gp++) {
            double fx = (nodeSign_[i][0] == nodeSign_[gp][0])? a : b;
            double fy = (nodeSign_[i][1] == nodeSign_[gp][1])? a : b;
            double fz = (nodeSign_[i][2] == nodeSign_[gp][2])? a : b;
            E(i, gp) = fx * fy * fz;
        }
    }
    return E;
}
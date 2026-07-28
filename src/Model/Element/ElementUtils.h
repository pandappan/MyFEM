// 构造单元在积分点处的数据：
// 1.形函数，形函数的全局导数，雅可比矩阵*权重
// 2.面，线外法向
#pragma once
#include "DenseMatrix.h"
#include "Quadrature.h"
#include "ShapeFun.h"


struct ElementIntegralData {
    unsigned int nen; // number of element nodes
    unsigned int ngs; // number of integrate points
    unsigned int nsd; // number of space dimensions

    // 形函数及其导数矩阵
    std::vector<std::vector<double>> N; // [numPoints] 每个积分点的形函数(numNodes)
    std::vector<DenseMatrix<double>> dN_dx; // [numPoints]全局坐标下的形函数的导数(dimCoord, numNode)

    std::vector<double> detJ_t_weight; // 积分权重
};

struct Normal {
    unsigned int numPoint;
    unsigned int dimCoord;
    DenseMatrix<double> normal; // (numPoint, dimCoord)
};

namespace ElementUtils {
    // 接受：形函数方案，积分方案，节点坐标数组，一次完成所有积分点运算
    // 输出：单元的Inform
    template<class ShapeFunType, class QuadFun>
    bool ComputeIntegralData(const DenseMatrix<double>& nodeCoords, ElementIntegralData& data);
    // 计算：(dx_i / dxi_j)^-1, detJ，一次计算一个积分点即可
    // 输入：单元节点坐标（dimCoord,numNode），形函数的局部导数，
    // 输出：单元雅可比逆矩阵以及行列式
    inline bool ComputeJacobian(const DenseMatrix<double>& nodeCoords, const DenseMatrix<double>& dN_dxi,
        DenseMatrix<double>& jacobianInverse, double& detJacobian);
    // 作用：计算形函数的全局导数，一次计算一个积分点数据
    // 输入：雅可比的逆矩阵，形函数的局部导数
    // 输出：dN_x
    inline bool ComputeGlobalDerivN(const DenseMatrix<double>& jacobianInverse,
        DenseMatrix<double>& dN_dxi, DenseMatrix<double>& dN_dx);
    // 作用：计算力学单元的B矩阵（1D，2D，3D）
    // 输入：形函数的全局导数
    // 输出：B矩阵(B1,B2,B3)排列形式
    inline bool MechanicsSolidBMatrix(const DenseMatrix<double>& dN_dx, DenseMatrix<double>& BMatrix);
    // 轴对称力学单元的B矩阵
    // 作用：计算轴对称单元的B矩阵
    // 输入：形函数的全局导数，形函数，半径
    inline bool MechanicsAxisBMatrix(const DenseMatrix<double>& dN_dx, const std::vector<double>& N,
        double radius, DenseMatrix<double>& Bmatrix);
    // 作用：计算线、面法向（兼容2D，3D连续介质单元计算的需求）
    // 输入：节点坐标（原始维度的节点坐标），形函数的局部导数
    // 输出：外法向
    inline bool SolidNormalDir(const DenseMatrix<double>& nodeCoords, const DenseMatrix<double>& dN_dxi,
        std::vector<double>& normVec);
}

template<class ShapeFunType, class QuadFun>
bool ElementUtils::ComputeIntegralData(const DenseMatrix<double>& nodeCoords, ElementIntegralData& data) {
    // 基础维度数据赋值，积分点数据
    data.nen = nodeCoords.GetCols();
    data.nsd = nodeCoords.GetRows();
    GaussPoints gps = QuadFun();
    data.ngs = gps.numPoints;

    // 设定结构体空间大小
    data.N.resize(data.ngs);
    data.dN_dx.resize(data.ngs);
    data.detJ_t_weight.resize(data.ngs);

    // 填充结构体数据
    for (unsigned int ip = 0; ip < data.ngs; ip++) {
        const auto& cha = gps.cha.Row(ip);
        const auto& weight = gps.weights[ip];
        // 形函数及其局部导数计算
        DenseMatrix<double> dN_dxi(data.nsd, data.nen);
        ShapeFun::ComputeShapFun<ShapeFunType>(cha, data.N[ip], dN_dxi);
        // 雅可比矩阵部分
        DenseMatrix<double> jacobianInverse(data.nsd, data.nsd);
        double detJacobian = 0.0;
        ComputeJacobian(nodeCoords, dN_dxi, jacobianInverse, detJacobian);
        // 全局导数
        ComputeGlobalDerivN(jacobianInverse, dN_dxi, data.dN_dx[ip]);
        // 雅可比行列式乘权重
        data.detJ_t_weight[ip] = detJacobian * weight;
    }
    return true;
}

inline bool ElementUtils::ComputeJacobian(const DenseMatrix<double> &nodeCoords,
    const DenseMatrix<double> &dN_dxi, DenseMatrix<double> &jacobianInverse, double& detJacobian) {
    const unsigned int nsd = nodeCoords.GetRows();
    const unsigned int nen = nodeCoords.GetCols();

    DenseMatrix<double> jacobian(nsd, nsd);
    jacobian.SetZero();
    for (unsigned int i = 0; i < nsd; i++) {
        for (unsigned int j = 0; j < nsd; j++) {
            double sum = 0.0;
            for (unsigned int I = 0; I < nen; I++) {
                sum += nodeCoords(i,I) * dN_dxi(j,I);
            }
            // 以转置形式设置雅可比矩阵
            jacobian(j,i) = sum;
        }
    }
    jacobianInverse = jacobian.Inverse();
    detJacobian = jacobian.Determinant();
    if (std::abs(detJacobian) < 1.0e-12) {
        throw std::runtime_error("ComputeJacobian: detJ is near zero!");
    }
    return true;
}

inline bool ElementUtils::ComputeGlobalDerivN(const DenseMatrix<double> &jacobianInverse,
    DenseMatrix<double> &dN_dxi, DenseMatrix<double> &dN_dx) {
    // 基础维度信息
    const unsigned int nsd = dN_dxi.GetRows();
    const unsigned int nen = dN_dxi.GetCols();

    // 初始化全局导数的空间
    dN_dx.Resize(nsd, nen);
    dN_dx.SetZero();
    // 计算全局导数
    for (unsigned int i = 0; i < nsd; i++) {
        for (unsigned int I = 0; I < nen; I++) {
            double sum = 0.0;
            for (unsigned int j = 0; j < nsd; j++) {
                sum += jacobianInverse(i,j) * dN_dxi(j, I);
            }
            dN_dx(i,I) = sum;
        }
    }
    return true;
}

inline bool ElementUtils::MechanicsSolidBMatrix(const DenseMatrix<double> &dN_dx, DenseMatrix<double> &BMatrix) {
    // 基础维度信息
    const unsigned int nsd = dN_dx.GetRows();
    const unsigned int nen = dN_dx.GetCols();
    const unsigned int nee = nsd * nen; // number of element equations

    // 初始化全局导数空间
    unsigned int nstrain = 0;
    switch (nsd) {
        case 1: nstrain = 1; break;
        case 2: nstrain = 3; break;
        case 3: nstrain = 6; break;
    }
    BMatrix.Resize(nstrain, nee);

    // 逐步计算节点B矩阵，填充整体B矩阵
    // 逻辑：根据局部方程号填充整体B矩阵
    if (nsd == 1) {
        for (unsigned int I = 0; I < nen; I++) {
            unsigned int colX = I * nsd;
            BMatrix(0,colX) = dN_dx(0, I);
        }
    } else if (nsd == 2) {
        for (unsigned int I = 0; I < nen; I++) {
            unsigned int colX = I * nsd;
            unsigned int colY = colX + 1;
            double dnx = dN_dx(0, I);
            double dny = dN_dx(1, I);
            BMatrix(0,colX) = dnx;
            BMatrix(1,colY) = dny;
            BMatrix(2,colX) = dny;
            BMatrix(2,colY) = dnx;
        }
    } else if (nsd == 3) {
        for (unsigned int I = 0; I < nen; I++) {
            unsigned int colX = I * nsd;
            unsigned int colY = colX + 1;
            unsigned int colZ = colY + 1;
            double dnx = dN_dx(0,I);
            double dny = dN_dx(1,I);
            double dnz = dN_dx(2,I);
            BMatrix(0,colX) = dnx;
            BMatrix(1,colY) = dny;
            BMatrix(2,colZ) = dnz;
            BMatrix(3,colX) = dny;
            BMatrix(3,colY) = dnx;
            BMatrix(4,colY) = dnz;
            BMatrix(4,colZ) = dny;
            BMatrix(5,colX) = dnz;
            BMatrix(5,colZ) = dnx;
        }
    }
    return true;
}

inline bool ElementUtils::MechanicsAxisBMatrix(const DenseMatrix<double>& dN_dx, const std::vector<double>& N,
    double radius, DenseMatrix<double>& BMatrix) {
    // 基础维度信息
    const unsigned int nsd = dN_dx.GetRows();
    const unsigned int nen = dN_dx.GetCols();
    const unsigned int nee = nsd * nen; // number of element equations
    BMatrix.Resize(4,nee);

    // 按照节点逐步填充B矩阵的值
    for (unsigned int I = 0; I < nen; I++) {
        unsigned int colR = I * 2;
        unsigned int colZ = colR + 1;
        const double dnr = dN_dx(0, I);
        const double dnz = dN_dx(1, I);
        const double dnt = N[I] / radius;
        BMatrix(0,colR) = dnr;
        BMatrix(1,colZ) = dnz;
        BMatrix(2,colR) = dnz;
        BMatrix(2,colZ) = dnr;
        BMatrix(3,colR) = dnt;
    }
    return true;
}

inline bool ElementUtils::SolidNormalDir(const DenseMatrix<double> &nodeCoords,
    const DenseMatrix<double> &dN_dxi, std::vector<double> &normVec) {
    // 确定基础维度信息，判断是1D还是2D外法向
    const unsigned int nsd = nodeCoords.GetRows();
    const unsigned int nen = nodeCoords.GetCols();
    // 不同维度的外法向计算
    if (nsd == 1) {
        normVec.resize(2);
        std::vector<double> t(2,0.0);
        for (unsigned int I = 0; I < nen; I++) {
            t[0] += nodeCoords(0,I) * dN_dxi(0,I);
            t[1] += nodeCoords(1,I) * dN_dxi(0,I);
        }
        normVec = {t[1],-t[0]};
    } else {
        normVec.resize(3);
        std::vector<double> txi(3,0.0);
        std::vector<double> teta(3,0.0);
        for (unsigned int I = 0; I < nen; I++) {
            txi[0] += nodeCoords(0,I) * dN_dxi(0,I);
            txi[1] += nodeCoords(1,I) * dN_dxi(0,I);
            txi[2] += nodeCoords(2,I) * dN_dxi(0,I);

            teta[0] += nodeCoords(0,I) * dN_dxi(1,I);
            teta[1] += nodeCoords(1,I) * dN_dxi(1,I);
            teta[2] += nodeCoords(2,I) * dN_dxi(1,I);
        }
        normVec[0] = txi[1] * teta[2] - txi[2] * teta[1];
        normVec[1] = txi[2] * teta[0] - txi[0] * teta[2];
        normVec[2] = txi[0] * teta[1] - txi[1] * teta[0];
    }
    // 外法向归一化
    double normVal = 0.0;
    for (const auto& e:normVec) {
        normVal += e*e;
    }
    normVal = std::sqrt(normVal);
    for (auto& e:normVec) {
        e /= normVal;
    }
    return true;
}
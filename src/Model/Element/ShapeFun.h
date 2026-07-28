#pragma once
#include <cassert>
#include <vector>
#include "../../Core/DenseMatrix.h"

// 类似FEAP的架构，将形函数及其全局导数作为一种公用方法，积分方案也作为一种公用数据
// 数值来源：FEAP shp1d.f / shp2d.f / shp3d.f / shp1dn.f / Shape.f90
namespace ShapeFun {
    // 两节点线元
    // 编号 1-2
    namespace L2 {
        constexpr unsigned int dimCoord = 1;
        constexpr unsigned int numNode = 2;
        inline void ComputeN(const std::vector<double> &cha, std::vector<double> &N) {
            assert(cha.size() == dimCoord);
            N.resize(numNode);
            const double xi = cha[0];
            N[0] = 0.5 * (1.0 - xi);
            N[1] = 0.5 * (1.0 + xi);
        }
        // dN_I/dxi_i，先维度再节点
        inline void ComputeDerivN(const std::vector<double> &cha, DenseMatrix<double> &dN_dxi) {
            assert(cha.size() == dimCoord);
            dN_dxi.Resize(dimCoord, numNode);
            dN_dxi(0,0) = -0.5;
            dN_dxi(0,1) =  0.5;
        }
    }
    // 三节点线元
    // 编号 1-3-2（节点3为中点，xi = -1, +1, 0）
    namespace L3 {
        constexpr unsigned int dimCoord = 1;
        constexpr unsigned int numNode = 3;
        inline void ComputeN(const std::vector<double> &cha, std::vector<double> &N) {
            assert(cha.size() == dimCoord);
            N.resize(numNode);
            const double xi = cha[0];
            N[0] = 0.5 * xi * (xi - 1.0);
            N[1] = 0.5 * xi * (xi + 1.0);
            N[2] = 1.0 - xi * xi;
        }
        // dN_I/dxi_i，先维度再节点
        inline void ComputeDerivN(const std::vector<double> &cha, DenseMatrix<double> &dN_dxi) {
            assert(cha.size() == dimCoord);
            dN_dxi.Resize(dimCoord, numNode);
            const double xi = cha[0];
            dN_dxi(0,0) = xi - 0.5;
            dN_dxi(0,1) = xi + 0.5;
            dN_dxi(0,2) = -2.0 * xi;
        }
    }
    // 四节点四边形单元
    // 编号
    // 4--3
    // |  |
    // 1--2
    namespace Q4 {
        constexpr unsigned int dimCoord = 2;
        constexpr unsigned int numNode = 4;
        inline void ComputeN(const std::vector<double> &cha, std::vector<double> &N) {
            assert(cha.size() == dimCoord);
            N.resize(numNode);
            const double xi = cha[0], eta = cha[1];
            N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
            N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
            N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
            N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
        }
        // dN_I/dxi_i，先维度再节点
        inline void ComputeDerivN(const std::vector<double> &cha, DenseMatrix<double> &dN_dxi) {
            assert(cha.size() == dimCoord);
            dN_dxi.Resize(dimCoord, numNode);
            const double xi = cha[0], eta = cha[1];
            dN_dxi(0,0) = - 0.25 * (1.0 - eta);
            dN_dxi(0,1) =   0.25 * (1.0 - eta);
            dN_dxi(0,2) =   0.25 * (1.0 + eta);
            dN_dxi(0,3) = - 0.25 * (1.0 + eta);

            dN_dxi(1,0) = - 0.25 * (1.0 - xi);
            dN_dxi(1,1) = - 0.25 * (1.0 + xi);
            dN_dxi(1,2) =   0.25 * (1.0 + xi);
            dN_dxi(1,3) =   0.25 * (1.0 - xi);
        }
    }
    // 8节点四边形单元（serendipity）
    // 编号
    // 4--7--3
    // |     |
    // 8     6
    // |     |
    // 1--5--2
    namespace Q8 {
        constexpr unsigned int dimCoord = 2;
        constexpr unsigned int numNode = 8;
        inline void ComputeN(const std::vector<double> &cha, std::vector<double> &N) {
            assert(cha.size() == dimCoord);
            N.resize(numNode);
            const double xi = cha[0], eta = cha[1];
            const double xi2 = 1.0 - xi * xi, eta2 = 1.0 - eta * eta;
            // 边中点
            N[4] = 0.5 * xi2 * (1.0 - eta);
            N[5] = 0.5 * (1.0 + xi) * eta2;
            N[6] = 0.5 * xi2 * (1.0 + eta);
            N[7] = 0.5 * (1.0 - xi) * eta2;
            // 角点 = 双线性 - 相邻两边中点的一半（FEAP Shap2D 的分级构造）
            N[0] = 0.25 * (1.0 - xi) * (1.0 - eta) - 0.5 * (N[7] + N[4]);
            N[1] = 0.25 * (1.0 + xi) * (1.0 - eta) - 0.5 * (N[4] + N[5]);
            N[2] = 0.25 * (1.0 + xi) * (1.0 + eta) - 0.5 * (N[5] + N[6]);
            N[3] = 0.25 * (1.0 - xi) * (1.0 + eta) - 0.5 * (N[6] + N[7]);
        }
        // dN_I/dxi_i，先维度再节点
        inline void ComputeDerivN(const std::vector<double> &cha, DenseMatrix<double> &dN_dxi) {
            assert(cha.size() == dimCoord);
            dN_dxi.Resize(dimCoord, numNode);
            const double xi = cha[0], eta = cha[1];
            const double xi2 = 1.0 - xi * xi, eta2 = 1.0 - eta * eta;
            // 边中点
            dN_dxi(0,4) = -xi * (1.0 - eta);
            dN_dxi(1,4) = -0.5 * xi2;
            dN_dxi(0,5) =  0.5 * eta2;
            dN_dxi(1,5) = -eta * (1.0 + xi);
            dN_dxi(0,6) = -xi * (1.0 + eta);
            dN_dxi(1,6) =  0.5 * xi2;
            dN_dxi(0,7) = -0.5 * eta2;
            dN_dxi(1,7) = -eta * (1.0 - xi);
            // 角点
            dN_dxi(0,0) = -0.25 * (1.0 - eta) - 0.5 * (dN_dxi(0,7) + dN_dxi(0,4));
            dN_dxi(1,0) = -0.25 * (1.0 - xi)  - 0.5 * (dN_dxi(1,7) + dN_dxi(1,4));
            dN_dxi(0,1) =  0.25 * (1.0 - eta) - 0.5 * (dN_dxi(0,4) + dN_dxi(0,5));
            dN_dxi(1,1) = -0.25 * (1.0 + xi)  - 0.5 * (dN_dxi(1,4) + dN_dxi(1,5));
            dN_dxi(0,2) =  0.25 * (1.0 + eta) - 0.5 * (dN_dxi(0,5) + dN_dxi(0,6));
            dN_dxi(1,2) =  0.25 * (1.0 + xi)  - 0.5 * (dN_dxi(1,5) + dN_dxi(1,6));
            dN_dxi(0,3) = -0.25 * (1.0 + eta) - 0.5 * (dN_dxi(0,6) + dN_dxi(0,7));
            dN_dxi(1,3) =  0.25 * (1.0 - xi)  - 0.5 * (dN_dxi(1,6) + dN_dxi(1,7));
        }
    }
    // 三节点三角形单元
    // 编号
    // 3
    // | \
    // 1--2
    namespace T3 {
        constexpr unsigned int dimCoord = 2;
        constexpr unsigned int numNode = 3;
        inline void ComputeN(const std::vector<double> &cha, std::vector<double> &N) {
            assert(cha.size() == dimCoord);
            N.resize(numNode);
            const double xi = cha[0], eta = cha[1], zeta = 1.0 - xi - eta;
            N[0] = zeta;
            N[1] = xi;
            N[2] = eta;
        }
        // dN_I/dxi_i，先维度再节点
        inline void ComputeDerivN(const std::vector<double> &cha, DenseMatrix<double> &dN_dxi) {
            assert(cha.size() == dimCoord);
            dN_dxi.Resize(dimCoord, numNode);
            // 注意：N = [1-xi-eta, xi, eta]
            dN_dxi(0,0) = -1.0;
            dN_dxi(0,1) =  1.0;
            dN_dxi(0,2) =  0.0;

            dN_dxi(1,0) = -1.0;
            dN_dxi(1,1) =  0.0;
            dN_dxi(1,2) =  1.0;
        }
    }
    // 6节点三角形单元
    // 编号（4:1-2中点, 5:2-3中点, 6:3-1中点）
    // 3
    // | \
    // 6   5
    // |    \
    // 1--4--2
    namespace T6 {
        constexpr unsigned int dimCoord = 2;
        constexpr unsigned int numNode = 6;
        inline void ComputeN(const std::vector<double> &cha, std::vector<double> &N) {
            assert(cha.size() == dimCoord);
            N.resize(numNode);
            const double L2 = cha[0], L3 = cha[1], L1 = 1.0 - L2 - L3;
            N[0] = L1 * (2.0 * L1 - 1.0);
            N[1] = L2 * (2.0 * L2 - 1.0);
            N[2] = L3 * (2.0 * L3 - 1.0);
            N[3] = 4.0 * L1 * L2;
            N[4] = 4.0 * L2 * L3;
            N[5] = 4.0 * L3 * L1;
        }
        // dN_I/dxi_i，先维度再节点
        inline void ComputeDerivN(const std::vector<double> &cha, DenseMatrix<double> &dN_dxi) {
            assert(cha.size() == dimCoord);
            dN_dxi.Resize(dimCoord, numNode);
            const double L2 = cha[0], L3 = cha[1], L1 = 1.0 - L2 - L3;
            // dL1/dxi = dL1/deta = -1
            dN_dxi(0,0) = -(4.0 * L1 - 1.0);
            dN_dxi(0,1) =   4.0 * L2 - 1.0;
            dN_dxi(0,2) =   0.0;
            dN_dxi(0,3) =   4.0 * (L1 - L2);
            dN_dxi(0,4) =   4.0 * L3;
            dN_dxi(0,5) = - 4.0 * L3;

            dN_dxi(1,0) = -(4.0 * L1 - 1.0);
            dN_dxi(1,1) =   0.0;
            dN_dxi(1,2) =   4.0 * L3 - 1.0;
            dN_dxi(1,3) = - 4.0 * L2;
            dN_dxi(1,4) =   4.0 * L2;
            dN_dxi(1,5) =   4.0 * (L1 - L3);
        }
    }
    // 四节点四面体单元
    namespace Tet4 {
        constexpr unsigned int dimCoord = 3;
        constexpr unsigned int numNode = 4;
        inline void ComputeN(const std::vector<double> &cha, std::vector<double> &N) {
            assert(cha.size() == dimCoord);
            N.resize(numNode);
            const double xi = cha[0], eta = cha[1], zeta = cha[2], omega = 1.0 - xi - eta - zeta;
            N[0] = omega;
            N[1] = xi;
            N[2] = eta;
            N[3] = zeta;
        }
        // dN_I/dxi_i，先维度再节点
        inline void ComputeDerivN(const std::vector<double> &cha, DenseMatrix<double> &dN_dxi) {
            assert(cha.size() == dimCoord);
            dN_dxi.Resize(dimCoord, numNode);
            dN_dxi(0,0) = -1.0;
            dN_dxi(0,1) =  1.0;
            dN_dxi(0,2) =  0.0;
            dN_dxi(0,3) =  0.0;

            dN_dxi(1,0) = -1.0;
            dN_dxi(1,1) =  0.0;
            dN_dxi(1,2) =  1.0;
            dN_dxi(1,3) =  0.0;

            dN_dxi(2,0) = -1.0;
            dN_dxi(2,1) =  0.0;
            dN_dxi(2,2) =  0.0;
            dN_dxi(2,3) =  1.0;
        }
    }
    // 10节点四面体单元
    // 编号：1-4为角点（同Tet4），棱中点 5:1-2, 6:2-3, 7:3-1, 8:1-4, 9:2-4, 10:3-4
    namespace Tet10 {
        constexpr unsigned int dimCoord = 3;
        constexpr unsigned int numNode = 10;
        inline void ComputeN(const std::vector<double> &cha, std::vector<double> &N) {
            assert(cha.size() == dimCoord);
            N.resize(numNode);
            const double L2 = cha[0], L3 = cha[1], L4 = cha[2], L1 = 1.0 - L2 - L3 - L4;
            N[0] = L1 * (2.0 * L1 - 1.0);
            N[1] = L2 * (2.0 * L2 - 1.0);
            N[2] = L3 * (2.0 * L3 - 1.0);
            N[3] = L4 * (2.0 * L4 - 1.0);
            N[4] = 4.0 * L1 * L2;
            N[5] = 4.0 * L2 * L3;
            N[6] = 4.0 * L3 * L1;
            N[7] = 4.0 * L1 * L4;
            N[8] = 4.0 * L2 * L4;
            N[9] = 4.0 * L3 * L4;
        }
        // dN_I/dxi_i，先维度再节点
        inline void ComputeDerivN(const std::vector<double> &cha, DenseMatrix<double> &dN_dxi) {
            assert(cha.size() == dimCoord);
            dN_dxi.Resize(dimCoord, numNode);
            const double L2 = cha[0], L3 = cha[1], L4 = cha[2], L1 = 1.0 - L2 - L3 - L4;
            const double g1 = 4.0 * L1 - 1.0;
            // d/dxi (xi = L2)
            dN_dxi(0,0) = -g1;
            dN_dxi(0,1) =  4.0 * L2 - 1.0;
            dN_dxi(0,2) =  0.0;
            dN_dxi(0,3) =  0.0;
            dN_dxi(0,4) =  4.0 * (L1 - L2);
            dN_dxi(0,5) =  4.0 * L3;
            dN_dxi(0,6) = -4.0 * L3;
            dN_dxi(0,7) = -4.0 * L4;
            dN_dxi(0,8) =  4.0 * L4;
            dN_dxi(0,9) =  0.0;
            // d/deta (eta = L3)
            dN_dxi(1,0) = -g1;
            dN_dxi(1,1) =  0.0;
            dN_dxi(1,2) =  4.0 * L3 - 1.0;
            dN_dxi(1,3) =  0.0;
            dN_dxi(1,4) = -4.0 * L2;
            dN_dxi(1,5) =  4.0 * L2;
            dN_dxi(1,6) =  4.0 * (L1 - L3);
            dN_dxi(1,7) = -4.0 * L4;
            dN_dxi(1,8) =  0.0;
            dN_dxi(1,9) =  4.0 * L4;
            // d/dzeta (zeta = L4)
            dN_dxi(2,0) = -g1;
            dN_dxi(2,1) =  0.0;
            dN_dxi(2,2) =  0.0;
            dN_dxi(2,3) =  4.0 * L4 - 1.0;
            dN_dxi(2,4) = -4.0 * L2;
            dN_dxi(2,5) =  0.0;
            dN_dxi(2,6) = -4.0 * L3;
            dN_dxi(2,7) =  4.0 * (L1 - L4);
            dN_dxi(2,8) =  4.0 * L2;
            dN_dxi(2,9) =  4.0 * L3;
        }
    }
    // 8节点六面体单元
    namespace H8 {
        constexpr unsigned int dimCoord = 3;
        constexpr unsigned int numNode = 8;
        inline void ComputeN(const std::vector<double> &cha, std::vector<double> &N) {
            assert(cha.size() == dimCoord);
            N.resize(numNode);
            double sxi[2]   = {1.0 - cha[0], 1.0 + cha[0]};
            double seta[2]  = {1.0 - cha[1], 1.0 + cha[1]};
            double szeta[2] = {1.0 - cha[2], 1.0 + cha[2]};
            static const int idx[8][3] = {
                {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
                {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
            };
            for (unsigned int I = 0; I < numNode; I++)
                N[I] = 0.125 * sxi[idx[I][0]] * seta[idx[I][1]] * szeta[idx[I][2]];
        }
        // dN_I/dxi_i，先维度再节点
        inline void ComputeDerivN(const std::vector<double> &cha, DenseMatrix<double> &dN_dxi) {
            assert(cha.size() == dimCoord);
            dN_dxi.Resize(dimCoord, numNode);
            double sxi[2]   = {1.0 - cha[0], 1.0 + cha[0]},   dxi[2] = {-1.0, 1.0};
            double seta[2]  = {1.0 - cha[1], 1.0 + cha[1]},  deta[2] = {-1.0, 1.0};
            double szeta[2] = {1.0 - cha[2], 1.0 + cha[2]}, dzeta[2] = {-1.0, 1.0};
            static const int idx[8][3] = {
                {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
                {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
            };
            for (unsigned int I = 0; I < numNode; I++) {
                dN_dxi(0, I) = 0.125 * dxi[idx[I][0]] * seta[idx[I][1]] * szeta[idx[I][2]];
                dN_dxi(1, I) = 0.125 * sxi[idx[I][0]] * deta[idx[I][1]] * szeta[idx[I][2]];
                dN_dxi(2, I) = 0.125 * sxi[idx[I][0]] * seta[idx[I][1]] * dzeta[idx[I][2]];
            }
        }
    }
    // 20节点六面体单元（serendipity）
    // 编号：1-8角点（同H8）；9-12底面棱中点(1-2,2-3,3-4,4-1)；
    //       13-16顶面棱中点(5-6,6-7,7-8,8-5)；17-20竖向棱中点(1-5,2-6,3-7,4-8)
    namespace H20 {
        constexpr unsigned int dimCoord = 3;
        constexpr unsigned int numNode = 20;
        // 各节点自然坐标
        inline const double (&NodeCoord())[20][3] {
            static const double xn[20][3] = {
                {-1,-1,-1}, { 1,-1,-1}, { 1, 1,-1}, {-1, 1,-1},
                {-1,-1, 1}, { 1,-1, 1}, { 1, 1, 1}, {-1, 1, 1},
                { 0,-1,-1}, { 1, 0,-1}, { 0, 1,-1}, {-1, 0,-1},
                { 0,-1, 1}, { 1, 0, 1}, { 0, 1, 1}, {-1, 0, 1},
                {-1,-1, 0}, { 1,-1, 0}, { 1, 1, 0}, {-1, 1, 0}
            };
            return xn;
        }
        inline void ComputeN(const std::vector<double> &cha, std::vector<double> &N) {
            assert(cha.size() == dimCoord);
            N.resize(numNode);
            const double xi = cha[0], eta = cha[1], zeta = cha[2];
            const auto &xn = NodeCoord();
            for (unsigned int I = 0; I < numNode; I++) {
                const double xI = xn[I][0], yI = xn[I][1], zI = xn[I][2];
                const double A = 1.0 + xi * xI, B = 1.0 + eta * yI, C = 1.0 + zeta * zI;
                if (I < 8)          // 角点
                    N[I] = 0.125 * A * B * C * (xi * xI + eta * yI + zeta * zI - 2.0);
                else if (xI == 0.0) // xi 方向棱中点
                    N[I] = 0.25 * (1.0 - xi * xi) * B * C;
                else if (yI == 0.0) // eta 方向棱中点
                    N[I] = 0.25 * A * (1.0 - eta * eta) * C;
                else                // zeta 方向棱中点
                    N[I] = 0.25 * A * B * (1.0 - zeta * zeta);
            }
        }
        // dN_I/dxi_i，先维度再节点
        inline void ComputeDerivN(const std::vector<double> &cha, DenseMatrix<double> &dN_dxi) {
            assert(cha.size() == dimCoord);
            dN_dxi.Resize(dimCoord, numNode);
            const double xi = cha[0], eta = cha[1], zeta = cha[2];
            const auto &xn = NodeCoord();
            for (unsigned int I = 0; I < numNode; I++) {
                const double xI = xn[I][0], yI = xn[I][1], zI = xn[I][2];
                const double A = 1.0 + xi * xI, B = 1.0 + eta * yI, C = 1.0 + zeta * zI;
                if (I < 8) {        // 角点
                    const double D = xi * xI + eta * yI + zeta * zI;
                    dN_dxi(0,I) = 0.125 * xI * B * C * (D + xi * xI - 1.0);
                    dN_dxi(1,I) = 0.125 * yI * A * C * (D + eta * yI - 1.0);
                    dN_dxi(2,I) = 0.125 * zI * A * B * (D + zeta * zI - 1.0);
                } else if (xI == 0.0) {
                    dN_dxi(0,I) = -0.5 * xi * B * C;
                    dN_dxi(1,I) =  0.25 * (1.0 - xi * xi) * yI * C;
                    dN_dxi(2,I) =  0.25 * (1.0 - xi * xi) * B * zI;
                } else if (yI == 0.0) {
                    dN_dxi(0,I) =  0.25 * xI * (1.0 - eta * eta) * C;
                    dN_dxi(1,I) = -0.5 * eta * A * C;
                    dN_dxi(2,I) =  0.25 * A * (1.0 - eta * eta) * zI;
                } else {
                    dN_dxi(0,I) =  0.25 * xI * B * (1.0 - zeta * zeta);
                    dN_dxi(1,I) =  0.25 * A * yI * (1.0 - zeta * zeta);
                    dN_dxi(2,I) = -0.5 * zeta * A * B;
                }
            }
        }
    }
    // 统一计算形函数及其局部导数
    template<class ShapeFunType>
    void ComputeShapFun(const std::vector<double> &cha, std::vector<double>& N, DenseMatrix<double>& dN_dxi) {
        ShapeFunType::ComputeN(cha,N);
        ShapeFunType::ComputeDerivN(cha, dN_dxi);
    }
}
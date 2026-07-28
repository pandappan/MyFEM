#pragma once
#include <cmath>
#include <vector>
#include "../../Core/DenseMatrix.h"

// 积分点数据
// cha:      (numPoints x dimCoord)，每行一个积分点的自然坐标
// weights:  积分权重
// 约定：
//   1) 四边形/六面体（Int1D/Int2D/Int3D）：标准 Gauss 权重，权和 = 2^dim
//   2) 三角形（TInt2D）：坐标为 (xi,eta)=(L2,L3)，权重已含面积因子 1/2，权和 = 1/2
//   3) 四面体（Tint3D）：坐标为 (xi,eta,zeta)=(L2,L3,L4)，权重已含体积因子 1/6，权和 = 1/6
//   故统一有  ∫f dΩ ≈ Σ w_i * f(cha_i) * detJ_i ，detJ 为 dN_dxi 算出的标准雅可比行列式
// 数值来源：FEAP int1d.f / int2d.f / int3d.f / tint2d.f / tint3d.f
struct GaussPoints {
    unsigned int dimCoord;
    unsigned int numPoints;
    DenseMatrix<double> cha;
    std::vector<double> weights;
};

namespace Quadrature {
    namespace Int1D {
        // 1维1点高斯积分
        inline GaussPoints Int1D1G() {
            GaussPoints data {};
            data.dimCoord = 1;
            data.numPoints = 1;
            data.cha.Resize(1,1);
            data.cha(0,0) = 0.0;
            data.weights.assign(1, 2.0);
            return data;
        }
        // 1维2点高斯积分
        inline GaussPoints Int1D2G() {
            GaussPoints data {};
            data.dimCoord = 1;
            data.numPoints = 2;
            data.cha.Resize(2,1);
            const double g = 1.0 / std::sqrt(3.0);
            data.cha(0,0) = -g;
            data.cha(1,0) =  g;
            data.weights.assign(2, 1.0);
            return data;
        }
        // 1维3点高斯积分
        inline GaussPoints Int1D3G() {
            GaussPoints data {};
            data.dimCoord = 1;
            data.numPoints = 3;
            data.cha.Resize(3,1);
            const double g = std::sqrt(0.6);
            data.cha(0,0) = -g;
            data.cha(1,0) =  0.0;
            data.cha(2,0) =  g;
            data.weights = { 5.0/9.0, 8.0/9.0, 5.0/9.0 };
            return data;
        }
    }
    namespace Int2D {
        // 2维四边形1点高斯积分
        inline GaussPoints Int2D1G() {
            GaussPoints data {};
            data.dimCoord = 2;
            data.numPoints = 1;
            data.cha.Resize(1,2);
            data.cha(0,0) = 0.0;
            data.cha(0,1) = 0.0;
            data.weights.assign(1, 4.0);
            return data;
        }
        // 2维四边形4点(2x2)高斯积分，点序与 FEAP int2d 一致
        inline GaussPoints Int2D4G() {
            GaussPoints data {};
            data.dimCoord = 2;
            data.numPoints = 4;
            data.cha.Resize(4,2);
            const double g = 1.0 / std::sqrt(3.0);
            static const int lr[4] = {-1, 1, 1,-1};
            static const int lz[4] = {-1,-1, 1, 1};
            for (unsigned int i = 0; i < 4; i++) {
                data.cha(i,0) = g * lr[i];
                data.cha(i,1) = g * lz[i];
            }
            data.weights.assign(4, 1.0);
            return data;
        }
        // 2维四边形9点(3x3)高斯积分，点序与 FEAP int2d 一致
        inline GaussPoints Int2D9G() {
            GaussPoints data {};
            data.dimCoord = 2;
            data.numPoints = 9;
            data.cha.Resize(9,2);
            const double g = std::sqrt(0.6);
            static const int lr[9] = {-1, 1, 1,-1, 0, 1, 0,-1, 0};
            static const int lz[9] = {-1,-1, 1, 1,-1, 0, 1, 0, 0};
            static const int lw[9] = {25,25,25,25,40,40,40,40,64};
            data.weights.resize(9);
            for (unsigned int i = 0; i < 9; i++) {
                data.cha(i,0) = g * lr[i];
                data.cha(i,1) = g * lz[i];
                data.weights[i] = lw[i] / 81.0;
            }
            return data;
        }
    }
    namespace TInt2D {
        // 2维三角形1点积分 O(h^2)
        inline GaussPoints TInt2D1G() {
            GaussPoints data {};
            data.dimCoord = 2;
            data.numPoints = 1;
            data.cha.Resize(1,2);
            const double third = 1.0/3.0;
            data.cha(0,0) = third;
            data.cha(0,1) = third;
            data.weights.assign(1, 0.5);
            return data;
        }
        // 2维三角形3点积分 O(h^3)，内点方案（FEAP tint2d l=-3）
        // 另一等精度方案为边中点：(0.5,0)、(0.5,0.5)、(0,0.5)，权同为 1/6
        inline GaussPoints TInt2D3G() {
            GaussPoints data {};
            data.dimCoord = 2;
            data.numPoints = 3;
            data.cha.Resize(3,2);
            const double a = 2.0/3.0, b = 1.0/6.0;
            // (xi,eta) = (L2,L3)
            data.cha(0,0) = b;  data.cha(0,1) = b;   // L1 = 2/3
            data.cha(1,0) = a;  data.cha(1,1) = b;   // L2 = 2/3
            data.cha(2,0) = b;  data.cha(2,1) = a;   // L3 = 2/3
            data.weights.assign(3, 1.0/6.0);
            return data;
        }
        // 2维三角形7点积分 O(h^5)
        inline GaussPoints TInt2D7G() {
            GaussPoints data {};
            data.dimCoord = 2;
            data.numPoints = 7;
            data.cha.Resize(7,2);
            data.weights.resize(7);
            const double third = 1.0/3.0;
            const double r0 = std::sqrt(15.0);
            const double r1 = 3.0/7.0;
            const double r2 = 2.0*r0/21.0;
            const double a = r1 + r2, b = 0.5*(1.0 - a);   // a≈0.797427, b≈0.101287
            const double c = r1 - r2, d = 0.5*(1.0 - c);   // c≈0.059716, d≈0.470142
            const double wa = 0.5 * (155.0 - r0)/1200.0;
            const double wc = 0.5 * (155.0 + r0)/1200.0;
            // 中心点
            data.cha(0,0) = third; data.cha(0,1) = third; data.weights[0] = 0.5*0.225;
            // L1=a / L2=a / L3=a
            data.cha(1,0) = b; data.cha(1,1) = b; data.weights[1] = wa;
            data.cha(2,0) = a; data.cha(2,1) = b; data.weights[2] = wa;
            data.cha(3,0) = b; data.cha(3,1) = a; data.weights[3] = wa;
            // L1=c / L2=c / L3=c
            data.cha(4,0) = d; data.cha(4,1) = d; data.weights[4] = wc;
            data.cha(5,0) = c; data.cha(5,1) = d; data.weights[5] = wc;
            data.cha(6,0) = d; data.cha(6,1) = c; data.weights[6] = wc;
            return data;
        }
    }
    namespace Int3D {
        // 三维六面体1点高斯积分
        inline GaussPoints Int3D1G() {
            GaussPoints data {};
            data.dimCoord = 3;
            data.numPoints = 1;
            data.cha.Resize(1,3);
            data.cha(0,0) = 0.0;
            data.cha(0,1) = 0.0;
            data.cha(0,2) = 0.0;
            data.weights.assign(1, 8.0);
            return data;
        }
        // 三维六面体8点(2x2x2)高斯积分
        inline GaussPoints Int3D8G() {
            GaussPoints data {};
            data.dimCoord = 3;
            data.numPoints = 8;
            data.cha.Resize(8,3);
            const double g = 1.0 / std::sqrt(3.0);
            static const int lr[4] = {-1, 1, 1,-1};
            static const int lz[4] = {-1,-1, 1, 1};
            for (unsigned int i = 0; i < 4; i++) {
                data.cha(i,0)   = g * lr[i];  data.cha(i,1)   = g * lz[i];  data.cha(i,2)   = -g;
                data.cha(i+4,0) = g * lr[i];  data.cha(i+4,1) = g * lz[i];  data.cha(i+4,2) =  g;
            }
            data.weights.assign(8, 1.0);
            return data;
        }
        // 三维六面体27点(3x3x3)高斯积分
        inline GaussPoints Int3D27G() {
            GaussPoints data {};
            data.dimCoord = 3;
            data.numPoints = 27;
            data.cha.Resize(27,3);
            data.weights.resize(27);
            const double g = std::sqrt(0.6);
            static const double s[3] = {-1.0, 0.0, 1.0};
            static const double w[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
            unsigned int p = 0;
            for (unsigned int k = 0; k < 3; k++)
                for (unsigned int j = 0; j < 3; j++)
                    for (unsigned int i = 0; i < 3; i++) {
                        data.cha(p,0) = g * s[i];
                        data.cha(p,1) = g * s[j];
                        data.cha(p,2) = g * s[k];
                        data.weights[p] = w[i] * w[j] * w[k];
                        p++;
                    }
            return data;
        }
    }
    namespace Tint3D {
        // 三维四面体1点积分 O(h^2)
        inline GaussPoints Tint3D1G() {
            GaussPoints data {};
            data.dimCoord = 3;
            data.numPoints = 1;
            data.cha.Resize(1,3);
            data.cha(0,0) = 0.25;
            data.cha(0,1) = 0.25;
            data.cha(0,2) = 0.25;
            data.weights.assign(1, 1.0/6.0);
            return data;
        }
        // 三维四面体4点积分 O(h^3)
        inline GaussPoints Tint3D4G() {
            GaussPoints data {};
            data.dimCoord = 3;
            data.numPoints = 4;
            data.cha.Resize(4,3);
            const double a = 0.5854101966249658;   // (5+3√5)/20
            const double b = 0.1381966011250105;   // (5-√5)/20
            // (xi,eta,zeta) = (L2,L3,L4)；第 i 点 L_i = a，其余 = b
            data.cha(0,0) = b; data.cha(0,1) = b; data.cha(0,2) = b;   // L1 = a
            data.cha(1,0) = a; data.cha(1,1) = b; data.cha(1,2) = b;
            data.cha(2,0) = b; data.cha(2,1) = a; data.cha(2,2) = b;
            data.cha(3,0) = b; data.cha(3,1) = b; data.cha(3,2) = a;
            data.weights.assign(4, 0.25/6.0);
            return data;
        }
        // 三维四面体11点积分 O(h^4)
        // 4 个角点 + 6 个棱中点 + 1 个中心点
        inline GaussPoints Tint3D11G() {
            GaussPoints data {};
            data.dimCoord = 3;
            data.numPoints = 11;
            data.cha.Resize(11,3);
            data.cha.SetZero();
            data.weights.resize(11);
            const double h = 0.5;
            // 角点 1-4：(L2,L3,L4)
            /* p0: 节点1 (0,0,0) 已置零 */
            data.cha(1,0) = 1.0;
            data.cha(2,1) = 1.0;
            data.cha(3,2) = 1.0;
            // 棱中点：1-2、2-3、3-4、1-4、2-4、1-3
            data.cha(4,0) = h;
            data.cha(5,0) = h; data.cha(5,1) = h;
            data.cha(6,1) = h; data.cha(6,2) = h;
            data.cha(7,2) = h;
            data.cha(8,0) = h; data.cha(8,2) = h;
            data.cha(9,1) = h;
            // 中心点
            data.cha(10,0) = 0.25; data.cha(10,1) = 0.25; data.cha(10,2) = 0.25;
            for (unsigned int i = 0; i < 4;  i++) data.weights[i] = 1.0/360.0;
            for (unsigned int i = 4; i < 10; i++) data.weights[i] = 1.0/90.0;
            data.weights[10] = 4.0/45.0;
            return data;
        }
    }
}
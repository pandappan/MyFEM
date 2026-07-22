//
// 轴对称 Q4 单元测试
//
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <vector>

#include "Model/Element/Q4_AX.h"
#include "Model/Element/CContinuumElement.h"
#include "Model/Node.h"
#include "Model/Material/CAxisymMaterial.h"
#include "Core/DenseMatrix.h"
#include "Core/Types.h"

namespace {

constexpr double PI = 3.14159265358979323846;

// 构造一个 [r1,r2] x [z1,z2] 的矩形轴对称单元
std::unique_ptr<CQ4_AX> MakeAxiElement(std::vector<CNode>& nodes,
                                       CAxisymMaterial& mat,
                                       double r1, double r2,
                                       double z1, double z2) {
    nodes.resize(4);
    // 局部节点顺序 N1(r1,z1) N2(r2,z1) N3(r2,z2) N4(r1,z2) 逆时针
    std::vector<double> c0{r1, z1, 0.0};
    std::vector<double> c1{r2, z1, 0.0};
    std::vector<double> c2{r2, z2, 0.0};
    std::vector<double> c3{r1, z2, 0.0};
    nodes[0].SetGeom(0, c0);
    nodes[1].SetGeom(1, c1);
    nodes[2].SetGeom(2, c2);
    nodes[3].SetGeom(3, c3);

    auto elem = std::unique_ptr<CQ4_AX>(new CQ4_AX());
    elem->SetElementType(ElementTypes::Q4_AX);

    std::vector<CNode*> nptr{&nodes[0], &nodes[1], &nodes[2], &nodes[3]};
    elem->SetupForTesting(nptr, &mat);
    elem->OnSetupComplete();   // 初始化积分点
    return elem;
}

} // namespace

// -------- 测试 1：材料弹性矩阵对称且正确 --------
TEST(Q4_AX, ElasticMatrix) {
    CAxisymMaterial mat;
    mat.E = 2.0e11;
    mat.nu = 0.3;

    DenseMatrix<double> D(4, 4);
    mat.ComputeElasticMatrix(D);

    double E = mat.E, nu = mat.nu;
    double a = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double b = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double G = E / (2.0 * (1.0 + nu));

    EXPECT_NEAR(D(0, 0), a, 1e-3);
    EXPECT_NEAR(D(0, 1), b, 1e-3);
    EXPECT_NEAR(D(2, 2), a, 1e-3);
    EXPECT_NEAR(D(3, 3), G, 1e-3);
    // 对称性
    for (unsigned int i = 0; i < 4; ++i)
        for (unsigned int j = 0; j < 4; ++j)
            EXPECT_NEAR(D(i, j), D(j, i), 1e-6);
}

// -------- 测试 2：单元体积应为 2πr̄·A --------
TEST(Q4_AX, VolumeIsRingVolume) {
    std::vector<CNode> nodes;
    CAxisymMaterial mat; mat.E = 1.0e7; mat.nu = 0.3;

    // r in [1,2], z in [0,1]，截面面积 A = 1
    auto elem = MakeAxiElement(nodes, mat, 1.0, 2.0, 0.0, 1.0);

    // 环体积 = 2π * r̄ * A，r̄ = 1.5，A = 1
    double expected = 2.0 * PI * 1.5 * 1.0;
    EXPECT_NEAR(elem->GetVolume(), expected, 1e-6);
}

// -------- 测试 3：刚度矩阵基本性质（对称、半正定、刚体位移零能量）--------
TEST(Q4_AX, StiffnessSymmetryAndRigidBody) {
    std::vector<CNode> nodes;
    CAxisymMaterial mat; mat.E = 2.0e11; mat.nu = 0.3;
    auto elem = MakeAxiElement(nodes, mat, 1.0, 2.0, 0.0, 1.0);

    DenseMatrix<double> Ke(8, 8);
    elem->ElementStiffness(Ke);

    // 对称性
    for (unsigned int i = 0; i < 8; ++i)
        for (unsigned int j = 0; j < 8; ++j)
            EXPECT_NEAR(Ke(i, j), Ke(j, i), std::abs(Ke(i, j)) * 1e-8 + 1e-3);

    // 对角元素为正
    for (unsigned int i = 0; i < 8; ++i)
        EXPECT_GT(Ke(i, i), 0.0);

    // 轴向刚体平移 uz=常数：εrr=εzz=εθθ=γrz=0，应变能应为0
    std::vector<double> uz(8, 0.0);
    for (unsigned int i = 0; i < 4; ++i) uz[2 * i + 1] = 1.0; // 所有 z 位移=1
    double energy = 0.0;
    for (unsigned int i = 0; i < 8; ++i)
        for (unsigned int j = 0; j < 8; ++j)
            energy += 0.5 * uz[i] * Ke(i, j) * uz[j];
    EXPECT_NEAR(energy, 0.0, 1e-3);

    // 注意：径向刚体平移 ur=常数 会产生 εθθ=ur/r ≠ 0，
    // 轴对称问题中这不是刚体模式，不应期望零能量。
}

// -------- 测试 4：B 矩阵 εθθ 行的正确性 --------
TEST(Q4_AX, BMatrixHoopStrain) {
    std::vector<CNode> nodes;
    CAxisymMaterial mat; mat.E = 1.0e7; mat.nu = 0.3;
    auto elem = MakeAxiElement(nodes, mat, 1.0, 2.0, 0.0, 1.0);

    // 均匀径向位移 ur=1
    for (int i = 0; i < 4; ++i) {
        nodes[i].displacement[UX] = 1.0;
        nodes[i].displacement[UY] = 0.0;
    }

    for (unsigned int ip = 0; ip < 4; ++ip) {
        auto strain = elem->ComputeStrainForTesting(ip);
        double r_ip = elem->GetRadiusAtIntegrationPoint(ip);
        // εrr = 0 (ur均匀), εθθ = ur/r = 1/r_ip
        EXPECT_NEAR(strain[0], 0.0, 1e-9);          // εrr
        EXPECT_NEAR(strain[2], 1.0 / r_ip, 1e-9);   // εθθ
        EXPECT_NEAR(strain[3], 0.0, 1e-9);          // γrz
    }
}

// -------- 测试 5：厚壁圆筒内压解析解对比（Lamé 解）--------
// 内半径 a，外半径 b，内压 p，外压 0
// 径向位移 u(r) = (1+nu)/E * p*a^2/(b^2-a^2) * [ (1-2nu) r + b^2/r ]
// 用多个轴对称单元离散一段圆筒，比较外/内表面径向位移
TEST(Q4_AX, ThickCylinderInternalPressure) {
    const double a = 1.0, b = 2.0;   // 内外半径
    const double p = 1.0e6;          // 内压
    const double E = 2.0e11, nu = 0.3;
    const double H = 0.1;            // 轴向高度（一层单元）

    CAxisymMaterial mat; mat.E = E; mat.nu = nu;

    const int NR = 20;               // 径向单元数
    const int NN_R = NR + 1;

    // 生成节点 (2 层 z=0, z=H)，节点编号：先内到外(z=0)，再内到外(z=H)
    std::vector<CNode> nodes(2 * NN_R);
    auto rOf = & [<sup>1</sup>](int i) { return a + (b - a) * i / double(NR); };
    for (int i = 0; i < NN_R; ++i) {
        std::vector<double> c0{rOf(i), 0.0, 0.0};
        std::vector<double> c1{rOf(i), H,   0.0};
        nodes[i].SetGeom(i, c0);
        nodes[NN_R + i].SetGeom(NN_R + i, c1);
    }

    // 生成单元
    std::vector<std::unique_ptr<CQ4_AX>> elems;
    for (int i = 0; i < NR; ++i) {
        auto e = std::unique_ptr<CQ4_AX>(new CQ4_AX());
        e->SetElementType(ElementTypes::Q4_AX);
        // N1(i,z0) N2(i+1,z0) N3(i+1,z1) N4(i,z1)
        std::vector<CNode*> nptr{
            &nodes[i], &nodes[i + 1],
            &nodes[NN_R + i + 1], &nodes[NN_R + i]};
        e->SetupForTesting(nptr, &mat);
        e->OnSetupComplete();
        elems.push_back(std::move(e));
    }

    // 组装全局刚度 (自由度 = 2 * 节点数)，用稠密矩阵便于测试
    unsigned int ndof = 2 * nodes.size();
    DenseMatrix<double> Kg(ndof, ndof);
    Kg.SetZero();
    for (auto& e : elems) {
        DenseMatrix<double> Ke(8, 8);
        e->ElementStiffness(Ke);
        const auto& en = e->GetNodes();
        // 局部->全局映射
        int gdof[8];
        for (int n = 0; n < 4; ++n) {
            gdof[2 * n]     = 2 * en[n]->Index;
            gdof[2 * n + 1] = 2 * en[n]->Index + 1;
        }
        for (int I = 0; I < 8; ++I)
            for (int J = 0; J < 8; ++J)
                Kg(gdof[I], gdof[J]) += Ke(I, J);
    }

    // 右端载荷：内表面 r=a 施加内压。等效节点力 F_r = p * (2πa) * H / 2 分配到内表面两节点
    std::vector<double> F(ndof, 0.0);
    double lineForce = p * 2.0 * PI * a * H; // 内表面总径向力
    F[2 * 0]            += lineForce / 2.0;  // 节点0 (r=a, z=0)
    F[2 * (NN_R + 0)]   += lineForce / 2.0;  // 节点NN_R (r=a, z=H)

    // 约束：z=0 与 z=H 两层的 uz 固定（平面应变段，防止轴向刚体位移）
    // 固定所有节点 uz，只保留径向自由度求解（一维径向问题）
    std::vector<bool> fixed(ndof, false);
    for (unsigned int n = 0; n < nodes.size(); ++n)
        fixed[2 * n + 1] = true; // 固定所有 uz

    // 罚函数法施加约束
    double penalty = 1.0e30;
    for (unsigned int i = 0; i < ndof; ++i)
        if (fixed[i]) { Kg(i, i) += penalty; F[i] = 0.0; }

    // 用 DenseMatrix 求解 Kg u = F
    DenseMatrix<double> Kinv = Kg.Inverse();
    std::vector<double> u = Kinv.DotVec(F);

    // 解析解
    auto uExact = & [<sup>2</sup>](double r) {
        return (1.0 + nu) / E * p * a * a / (b * b - a * a)
               * ((1.0 - 2.0 * nu) * r + b * b / r);
    };

    // 比较内表面 (r=a, 节点0) 与外表面 (r=b, 节点NR) 的径向位移
    double ur_inner = u[2 * 0];
    double ur_outer = u[2 * NR];

    EXPECT_NEAR(ur_inner, uExact(a), uExact(a) * 0.03); // 3% 容差
    EXPECT_NEAR(ur_outer, uExact(b), uExact(b) * 0.03);
}
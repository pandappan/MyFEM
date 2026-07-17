#include <gtest/gtest.h>
#include "Node.h"
#include "Material/BarMaterial.h"
#include "Material/CPlaneStressMaterial.h"
#include "Element/Bar3D.h"
#include "Element/Q4.h"

TEST(BodyForc, Bar3DGravityHalvedBetweenNodes) {
    // 基本参数
    const double L = 2.0, Area = 1.0, rho = 100.0;
    // 节点参数
    std::vector<CNode> nodes;
    nodes.emplace_back(0.0,0.0,0.0);
    nodes.emplace_back(L,0.0,0.0);
    for (unsigned int i = 0; i < nodes.size(); i++) {
        nodes[i].Index = i;
    }
    // 材料参数
    auto mat = std::unique_ptr<CBarMaterial>(new CBarMaterial());
    mat->rho = rho;
    mat->Area = Area;
    // 单元参数
    auto elem = std::unique_ptr<CBar3D>(new CBar3D());
    std::vector<CNode*> nds = {&nodes[0],&nodes[1]};
    elem->SetupForTesting(nds,mat.get());
    // 体力加速度
    const double g[3] = {0.0, -10.0, 0.0};
    elem->CalculateBodyForce(g);
    // 验证
    EXPECT_NEAR(nodes[0].NodeForce[0],0.0,1.0e-8);
    EXPECT_NEAR(nodes[0].NodeForce[1],-1000.0,1.0e-8);
    EXPECT_NEAR(nodes[0].NodeForce[2],0.0,1.0e-8);
    EXPECT_NEAR(nodes[1].NodeForce[0],0.0,1.0e-8);
    EXPECT_NEAR(nodes[1].NodeForce[1],-1000.0,1.0e-8);
    EXPECT_NEAR(nodes[1].NodeForce[2],0.0,1.0e-8);
}

TEST(BodyForce, Q4UniformGravity) {
    // 单位方形 Q4，ρ=100, thk=1
    // 重力 g = [0, -10, 0]
    // 单元总重 = ρ·V·|g| = 100·1·10 = 1000
    // 应均匀分给 4 个节点（因为几何对称），每节点 -250 in Y

    std::vector<CNode> nodes;
    nodes.emplace_back(0.0, 0.0, 0.0);
    nodes.emplace_back(1.0, 0.0, 0.0);
    nodes.emplace_back(1.0, 1.0, 0.0);
    nodes.emplace_back(0.0, 1.0, 0.0);
    for (unsigned int i = 0; i < 4; ++i) nodes[i].Index = i;

    auto mat = std::unique_ptr<CPlaneStressMaterial>(new CPlaneStressMaterial());
    mat->E = 1e6; mat->nu = 0.3; mat->thk = 1.0; mat->rho = 100.0;

    auto elem = std::unique_ptr<CQ4>(new CQ4());
    std::vector<CNode*> nps = {&nodes[0], &nodes[1], &nodes[2], &nodes[3]};
    elem->SetupForTesting(nps, mat.get());
    elem->InitializeIntegrationPoints();

    double g[3] = {0.0, -10.0, 0.0};
    elem->CalculateBodyForce(g);

    // 4 个节点各分 -250 in Y
    // 局部 dof: [UX_N1, UY_N1, UX_N2, UY_N2, ...]
    double total_y = 0.0;
    for (unsigned int i = 0; i < 4; ++i) {
        EXPECT_NEAR(nodes[i].NodeForce[0],     0.0,   1e-10) << "Node " << i+1 << " UX";
        EXPECT_NEAR(nodes[i].NodeForce[1], -250.0, 1e-10) << "Node " << i+1 << " UY";
        total_y += nodes[i].NodeForce[1];
    }

    // 总力守恒
    EXPECT_NEAR(total_y, -1000.0, 1e-10);
}
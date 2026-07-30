#include <gtest/gtest.h>
#include "DenseMatrix.h"
#include "Node.h"
#include "Element/Q4.h"
#include "Material/PlaneStressMaterial.h"
// 使用一个梯形单元（畸变的 Q4）
// 施加沿 x 方向的纯应变场 u = ε·x, v = 0
// 期望：ε_xx = ε 常数，ε_yy = 0，γ_xy = 0
//
// 这个测试如果 Jacobian 转置错了会挂
TEST(Q4Distorted, PureStrainOnTrapezoid) {
    std::vector<Node> nodes;
    // 梯形：底边 [0,2]，顶边 [0.3, 1.7]
    nodes.emplace_back(0.0, 0.0, 0.0);
    nodes.emplace_back(2.0, 0.0, 0.0);
    nodes.emplace_back(1.7, 1.0, 0.0);
    nodes.emplace_back(0.3, 1.0, 0.0);
    for (unsigned int i = 0; i < 4; ++i) nodes[i].Index = i + 1;
    
    auto mat = std::unique_ptr<PlaneStressMaterial>(new PlaneStressMaterial());
    mat->E = 1.0; mat->nu = 0.0; mat->thk = 1.0;
    
    auto elem = std::unique_ptr<Q4>(new Q4());
    std::vector<Node*> nps = {&nodes[0], &nodes[1], &nodes[2], &nodes[3]};
    elem->SetupForTesting(nps, mat.get());
    elem->InitializeIntegrationPoints();
    
    // 手动施加纯 x 应变：u = 0.01·x, v = 0
    const double eps = 0.01;
    for (auto& n : nodes) {
        n.displacement[0] = eps * n.XYZ[0];
        n.displacement[1] = 0.0;
    }
    
    // 检查每个积分点的应变
    for (unsigned int ip = 0; ip < elem->GetNumIntegrationPoints(); ++ip) {
        auto strain = elem->ComputeStrainAtIntegrationPoint(ip);
        EXPECT_NEAR(strain[0], eps, 1e-10) << "eps_xx at ip " << ip;
        EXPECT_NEAR(strain[1], 0.0, 1e-10) << "eps_yy at ip " << ip;
        EXPECT_NEAR(strain[2], 0.0, 1e-10) << "gamma_xy at ip " << ip;
    }
}
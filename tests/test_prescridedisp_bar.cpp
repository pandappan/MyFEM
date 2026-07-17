#include <gtest/gtest.h>
#include "Model.h"
#include "Node.h"
#include "Element/Bar3D.h"
#include "Material/BarMaterial.h"
#include "Assembly.h"

TEST(IntegrationTest, BarWithPrescribedDisplacement) {
    // Bar: L = 1, EA = 100
    // node 1固定， nodes 强制UX=0.02
    // sigma = E * epsilon = E * (UX/L) = 2
    // R1 = -2
    // R2 = 2
    const double E = 100.0, A = 1.0, L = 1.0;
    const double u_prescribed = 0.02;
    Model model;
    model.dimension = 3;
    model.nodes.emplace_back(0.0, 0.0, 0.0);
    model.nodes.emplace_back(L, 0.0, 0.0);
    model.nodes[0].Index = 1;
    model.nodes[1].Index = 2;
    model.nodes[0].bcode[UX] = 1; // 全约束
    model.nodes[0].bcode[UY] = 1;
    model.nodes[0].bcode[UZ] = 1;
    model.nodes[1].bcode[UX] = 2; // 给定2指定位移，其余约束
    model.nodes[1].bcode[UY] = 1;
    model.nodes[1].bcode[UZ] = 1;
    model.nodes[1].SetPreDisp(UX, u_prescribed);
    // -------- 建 Element Group --------
    CElementGroup group;
    group.SetTypeForTesting(ElementTypes::Bar3D);

    // 材料
    auto mat = std::unique_ptr<CBarMaterial>(new CBarMaterial());
    mat->nset = 1;
    mat->E    = E;
    mat->Area = A;
    model.materials.push_back(std::move(mat));

    // 单元
    auto elem = std::unique_ptr<CBar3D>(new CBar3D());
    std::vector<CNode*> nodes = {&model.nodes[0], &model.nodes[1]};
    std::vector<unsigned int> connectivity = {0, 1};
    elem->SetElementInfo(0, model.materials[0].get(), connectivity, model.nodes);
    group.AddElement(std::move(elem));
    model.groups.push_back(std::move(group));

    Assembler::CalculateEquationNumber(model);
    ASSERT_EQ(model.neq, 0) << "All dofs are constrained, system has no free dof";
    Assembler::CalculateLocationMatrix(model);
    Assembler::CalculateNodalBCForce(model);

    // 验证反力
    const double expect_rf = E * A * u_prescribed / L;
    EXPECT_NEAR(nodes[0]->NodeBCForce[0], -expect_rf, 1e-8);
    EXPECT_NEAR(nodes[1]->NodeBCForce[0], +expect_rf, 1e-8);
}
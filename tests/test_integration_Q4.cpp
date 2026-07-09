//
// Created by Administrator on 2026/7/8.
//
#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include "Model.h"
#include "Node.h"
#include "Element/Q4.h"
#include "Element/ElementGroup.h"
#include "Material/CPlaneStressMaterial.h"
#include "Assembly.h"
#include "Solver.h"
TEST(IntegrationTest, Q4_UniaxialTension_UsingFullPipeline) {
    // 单个 Q4 单元，纯 X 向拉伸
    // 单位方形，E=1e6, ν=0, thk=1, F=1
    // 期望：u_x(N2)=u_x(N3)=F/E, u_y=0（ν=0 所以横向无位移）

    const double E = 1.0e6;
    const double nu = 0.0;
    const double F_per_node = 0.5;   // Node 2 和 3 各 0.5
    const double expected_ux = 2 * F_per_node / E;   // σ=1/1=1, ε=1/E

    Model model;
    model.dimension = 2;

    // 4 节点：单位方形
    model.nodes.emplace_back(0.0, 0.0, 0.0);
    model.nodes.emplace_back(1.0, 0.0, 0.0);
    model.nodes.emplace_back(1.0, 1.0, 0.0);
    model.nodes.emplace_back(0.0, 1.0, 0.0);
    for (unsigned int i = 0; i < 4; ++i)
        model.nodes[i].NodeNumber = i + 1;

    // 2D 分析：UZ 全部约束
    for (auto& n : model.nodes) n.bcode[UZ] = 1;

    // Node 1: 全固定
    model.nodes[0].bcode[UX] = 1;
    model.nodes[0].bcode[UY] = 1;
    // Node 2: X 自由，Y 固定
    model.nodes[1].bcode[UY] = 1;
    // Node 3: 全自由
    // Node 4: X 固定，Y 自由
    model.nodes[3].bcode[UX] = 1;

    // 载荷
    model.nodes[1].NodeForce[UX] = F_per_node;
    model.nodes[2].NodeForce[UX] = F_per_node;

    // Group
    CElementGroup group;
    group.SetTypeForTesting(ElementTypes::Q4_PS);

    auto mat = std::unique_ptr<CPlaneStressMaterial>(new CPlaneStressMaterial());
    mat->nset = 1;
    mat->E   = E;
    mat->nu  = nu;
    mat->thk = 1.0;
    group.AddMaterialForTesting(std::move(mat));

    auto elem = std::unique_ptr<CQ4>(new CQ4());
    std::vector<CNode*> nodes = {&model.nodes[0], &model.nodes[1], &model.nodes[2], &model.nodes[3]};
    elem->SetupForTesting(
        nodes,
        &group.GetMaterial(0)
    );
    elem.get()->InitializeIntegrationPoints();
    group.AddElementForTesting(std::move(elem));

    model.groups.push_back(std::move(group));

    // Pipeline
    Assembler::CalculateEquationNumber(model);
    Assembler::CalculateLocationMatrix(model);
    Assembler::AllocateLinearSystem(model);
    Assembler::AssembleForce(model);
    Assembler::AssembleStiffnessAndConstraintCorrection(model);

    CLDLTSolver solver(*model.K);
    solver.LDLT();
    solver.BackSubstitution(model.force);
    Assembler::WriteDisplacementToNodes(model);

    // 断言
    EXPECT_NEAR(model.nodes[1].Displacement[UX], expected_ux, 1e-10);
    EXPECT_NEAR(model.nodes[2].Displacement[UX], expected_ux, 1e-10);
    EXPECT_NEAR(model.nodes[2].Displacement[UY], 0.0, 1e-10);
    EXPECT_NEAR(model.nodes[3].Displacement[UY], 0.0, 1e-10);
}
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
TEST(StressSmoothing, Q4UniformStressField) {
    // 单个 Q4 单元，纯 X 向拉伸
    // 单位方形，E=1e6, ν=0, thk=1, F=1
    // 期望：u_x(N2)=u_x(N3)=F/E, u_y=0（ν=0 所以横向无位移）

    const double E = 1.0e6;
    const double nu = 0.0;
    const double F_per_node = 0.5;   // Node 2 和 3 各 0.5

    Model model;
    model.dimension = 2;

    // 4 节点：单位方形
    model.nodes.emplace_back(0.0, 0.0, 0.0);
    model.nodes.emplace_back(1.0, 0.0, 0.0);
    model.nodes.emplace_back(1.0, 1.0, 0.0);
    model.nodes.emplace_back(0.0, 1.0, 0.0);
    for (unsigned int i = 0; i < 4; ++i)
        model.nodes[i].Index = i;

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
    model.nodes[1].nodeForce[UX] = F_per_node;
    model.nodes[2].nodeForce[UX] = F_per_node;

    // Group
    CElementGroup group;
    group.SetTypeForTesting(ElementTypes::Q4_PS);

    auto mat = std::unique_ptr<CPlaneStressMaterial>(new CPlaneStressMaterial());
    mat->nset = 1;
    mat->E   = E;
    mat->nu  = nu;
    mat->thk = 1.0;
    model.materials.push_back(std::move(mat));

    auto elem = std::unique_ptr<CQ4>(new CQ4());
    std::vector<CNode*> nodes = {&model.nodes[0], &model.nodes[1], &model.nodes[2], &model.nodes[3]};
    std::vector<unsigned int> connectivity = {0,1,2,3};
    elem->SetElementInfo(0, model.materials[0].get(), connectivity, model.nodes);
    elem->InitializeIntegrationPoints();
    group.AddElement(std::move(elem));

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
    Assembler::CalculateNodalStress(model);


    // 均匀应力外推也是均匀
    for (unsigned int i = 0; i < 4; i++) {
        ASSERT_EQ(model.nodes[i].stress.size(),3);
        EXPECT_NEAR(model.nodes[i].stress[0], 1.0, 1e-8) << "Node " << i+1 << " sxx";
        EXPECT_NEAR(model.nodes[i].stress[1], 0.0, 1e-8) << "Node " << i+1 << " syy";
        EXPECT_NEAR(model.nodes[i].stress[2], 0.0, 1e-8) << "Node " << i+1 << " sxy";
    }
}


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

TEST(IntegrationTest, Q4SurfaceLoadTopEdge) {
    const double E = 100.0;
    const double q = 10.0;           // 顶边均布拉力
    const double expected_uy = q / E; // 拉伸位移 = q·L/(E·A)，L=1, A=1 → q/E

    Model model;
    model.dimension = 2;

    // 单位方形
    model.nodes.emplace_back(0.0, 0.0, 0.0);
    model.nodes.emplace_back(1.0, 0.0, 0.0);
    model.nodes.emplace_back(1.0, 1.0, 0.0);
    model.nodes.emplace_back(0.0, 1.0, 0.0);
    for (unsigned int i = 0; i < 4; ++i)
        model.nodes[i].Index = i;
    
    // 2D: UZ 全约束
    for (auto& n : model.nodes) n.bcode[UZ] = 1;
    
    // 约束：底边 UY=0，N1 的 UX 也=0（防止水平刚体移动）
    model.nodes[0].bcode[UX] = 1;
    model.nodes[0].bcode[UY] = 1;
    model.nodes[1].bcode[UY] = 1;
    // N3, N4 全部自由（除了 UZ）
    
    // Group + Material + Element
    CElementGroup group;
    group.SetTypeForTesting(ElementTypes::Q4_PS);
    
    auto mat = std::unique_ptr<CPlaneStressMaterial>(new CPlaneStressMaterial());
    mat->nset = 1;  mat->E = E;  mat->nu = 0.0;  mat->thk = 1.0;
    model.materials.push_back(std::move(mat));
    
    auto elem = std::unique_ptr<CQ4>(new CQ4());
    elem->SetElementNumber(1);   // ★ 面力查找需要
    std::vector<CNode*> nps = {&model.nodes[0], &model.nodes[1],
                               &model.nodes[2], &model.nodes[3]};
    std::vector<unsigned int> connectivity = {0,1,2,3};
    elem->SetElementInfo(0,model.GetMaterialPtr(0),connectivity,model.nodes);
    elem->InitializeIntegrationPoints();
    group.AddElement(std::move(elem));
    model.groups.push_back(std::move(group));
    
    //面力：顶边 (faceID=3, N3-N4, 1-based=3) 施加 UY 方向 q=10
    SurfaceLoad sl;
    sl.elemId_0 = 0;
    sl.faceId_0 = 2;    // 1-based，对应顶边
    sl.dof_0    = 1;
    sl.value  = q;
    model.sloads.push_back(sl);
    
    // Pipeline
    Assembler::CalculateEquationNumber(model);
    Assembler::CalculateLocationMatrix(model);
    Assembler::AllocateLinearSystem(model);
    Assembler::InitializeElementMap(model);
    Assembler::ConvertSLoadsToCLoads(model);
    Assembler::AssembleForce(model);
    Assembler::AssembleStiffnessAndConstraintCorrection(model);
    
    CLDLTSolver solver(*model.K);
    solver.LDLT();
    solver.BackSubstitution(model.force);
    Assembler::WriteDisplacementToNodes(model);
    
    // 断言：N3 和 N4 的 UY 都是 q/E
    EXPECT_NEAR(model.nodes[2].displacement[UY], expected_uy, 1e-8);
    EXPECT_NEAR(model.nodes[3].displacement[UY], expected_uy, 1e-8);
    
    // 因为 ν=0，UX 应该几乎为 0
    EXPECT_NEAR(model.nodes[2].displacement[UX], 0.0, 1e-8);
}
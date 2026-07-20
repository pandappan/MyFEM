#include <gtest/gtest.h>
#include <memory>
#include "Node.h"
#include "Element/Q4.h"
#include "Material/CPlaneStressMaterial.h"

TEST(NodeDofTest, DefaultMaskIsZero) {
    CNode node;
    EXPECT_EQ(node.activeMask, 0);
    EXPECT_FALSE(node.HasAnyActiveDof());
}

TEST(NodeDofTest, ActiveDof) {
    CNode node;
    node.ActivateDof(UX);
    node.ActivateDof(ROTZ);
    EXPECT_TRUE(node.IsDofActive(UX));
    EXPECT_TRUE(node.IsDofActive(ROTZ));
    EXPECT_FALSE(node.IsDofActive(UY));
    EXPECT_FALSE(node.IsDofActive(UZ));
}

TEST(NodeDofTest, GenerationEquation) {
    CNode node;
    node.ActivateDof(UX);
    node.ActivateDof(ROTZ);
    unsigned int neq = 0;
    node.GenerateNodeEquation(neq);
    EXPECT_EQ(neq, 2);
    EXPECT_EQ(node.eqn[0], 1);
    EXPECT_EQ(node.eqn[1], 0);
    EXPECT_EQ(node.eqn[2], 0);
    EXPECT_EQ(node.eqn[3], 0);
    EXPECT_EQ(node.eqn[4], 0);
    EXPECT_EQ(node.eqn[5], 2);
    EXPECT_EQ(node.eqn[6], 0);
}

TEST(NodeDofTest, ElementRegisterNodeDof) {
    std::vector<CNode> nodes(4);
    std::unique_ptr<CElement> elem = std::make_unique<CQ4>();
    std::unique_ptr<CMaterial> mat = std::make_unique<CPlaneStressMaterial>();

    std::vector<unsigned int> connectivity = {0,1,2,3};
    elem->SetElementInfo(0,mat.get(),connectivity,nodes);

    unsigned int neq = 0;
    for (auto& node : nodes) {
        node.GenerateNodeEquation(neq);
    }
    // 所有节点的UX，UY被激活,UZ关闭
    EXPECT_EQ(neq, 8);
    EXPECT_EQ(nodes[0].eqn[0], 1);
    EXPECT_EQ(nodes[0].eqn[1], 2);
    EXPECT_FALSE(nodes[0].IsDofActive(2));
}
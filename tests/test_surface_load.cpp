#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include "Node.h"
#include "Element/Q4.h"
#include "Material/PlaneStressMaterial.h"

class Q4SurfaceLoadFixture : public :: testing::Test {
protected:
    std::vector<Node> nodes_;
    std::unique_ptr<PlaneStressMaterial> mat_;
    std::unique_ptr<Q4> elem_;
    void SetUp() {
        nodes_.emplace_back(0.0,0.0,0.0);
        nodes_.emplace_back(2.0,0.0,0.0);
        nodes_.emplace_back(2.0,3.0,0.0);
        nodes_.emplace_back(0.0,3.0,0.0);

        mat_.reset(new PlaneStressMaterial());
        mat_->E = 1e6;
        mat_->nu = 0.3;
        mat_->thk = 1.0;

        elem_.reset(new Q4());
        std::vector<Node*> nps = {&nodes_[0], &nodes_[1], &nodes_[2], &nodes_[3]};
        elem_->SetupForTesting(nps,mat_.get());
    }
};

TEST_F(Q4SurfaceLoadFixture, Surface0UniformLoadX) {
    elem_->CalculateSurfaceLoad(0,UX,5);
    EXPECT_NEAR(nodes_[0].nodeForce[0],5.0,1e-10);
    EXPECT_NEAR(nodes_[0].nodeForce[1],0.0,1e-10);
    EXPECT_NEAR(nodes_[1].nodeForce[0],5.0,1e-10);
    EXPECT_NEAR(nodes_[1].nodeForce[1],0.0,1e-10);
}
// test_face_mapping.cpp
#include <gtest/gtest.h>
#include "Element/Q4.h"

TEST(Q4FaceMapping, AllFourFacesMapCorrectly) {
    // Directly call the mapping function (make it public or friend)
    Q4 q4;
    // 设置4个哑节点

    auto face0 = q4.GetFaceNodesLocalID(0);
    EXPECT_EQ(face0[0], 0);
    EXPECT_EQ(face0[1], 1);

    auto face1 = q4.GetFaceNodesLocalID(1);
    EXPECT_EQ(face1[0], 1);
    EXPECT_EQ(face1[1], 2);

    auto face2 = q4.GetFaceNodesLocalID(2);
    EXPECT_EQ(face2[0], 2);
    EXPECT_EQ(face2[1], 3);

    auto face3 = q4.GetFaceNodesLocalID(3);
    EXPECT_EQ(face3[0], 3);
    EXPECT_EQ(face3[1], 0) << "Face 3 should close the quad: N4->N1";
}

TEST(Q4FaceMapping, InvalidFaceThrows) {
    Q4 q4;
    EXPECT_THROW(q4.GetFaceNodesLocalID(4), std::out_of_range);
    EXPECT_THROW(q4.GetFaceNodesLocalID(99), std::out_of_range);
}

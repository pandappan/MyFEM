#include <gtest/gtest.h>
#include "../Core/DenseMatrix.h"

TEST(DenseMatrix, DefaultConstruct) {
    DenseMatrix<double> m;
    EXPECT_EQ(m.GetRows(), 0u);
    EXPECT_EQ(m.GetCols(), 0u);
}

TEST(DenseMatrix, Determinant2x2) {
    DenseMatrix<double> m(2,2);
    m(0,0) = 1.0; m(0,1) = 2.0;
    m(1,0) = 3.0; m(1,1) = 4.0;
    EXPECT_DOUBLE_EQ(m.Determinant(), -2.0);
}

TEST(DenseMatrix, Determinate3x3) {
    DenseMatrix<double> m(3,3);
    m(0,0) = 1.0; m(0,1) = 2.0; m(0,2) = 3.0;
    m(1,0) = 4.0; m(1,1) = 5.0; m(1,2) = 6.0;
    m(2,0) = 7.0; m(2,1) = 8.0; m(2,2) = 9.0;
    EXPECT_DOUBLE_EQ(m.Determinant(), 0.0);
}

TEST(DenseMatrix, InverseTimesOriginalIsIdentity) {
    DenseMatrix<double> m(4,4);
    DenseMatrix<double> Identity(4, 4);
    Identity.SetZero();
    double values[16] = {3, 8, 5, 6, 6, 9, 6, 2, 13, 8, 13, 16, 8, 13, 8, 8};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            m(i, j) = values[i*4 + j];
            if (i==j) {
                Identity(i,j) = 1.0;
            }
        }
    }
    DenseMatrix<double> Invm = m.Inverse();
    DenseMatrix<double> InverseTimesOriginal = m.DotMat(Invm);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double expected = (i==j) ? 1.0 : 0.0;
            EXPECT_NEAR(InverseTimesOriginal(i, j), expected, 1.0e-10);
        }
    }
}

TEST(DenseMatrix, TransposeTwiceIsIdentity) {
    DenseMatrix<double> m(2, 3);
    m(0,0)=1; m(0,1)=2; m(0,2)=3;
    m(1,0)=4; m(1,1)=5; m(1,2)=6;

    DenseMatrix<double> mt = m.Transpose();
    EXPECT_EQ(mt.GetRows(), 3u);
    EXPECT_EQ(mt.GetCols(), 2u);
    EXPECT_DOUBLE_EQ(mt(0,0), 1.0);
    EXPECT_DOUBLE_EQ(mt(1,0), 2.0);
    EXPECT_DOUBLE_EQ(mt(2,1), 6.0);

    // 转置两次回到原矩阵
    DenseMatrix<double> mtt = mt.Transpose();
    for (unsigned int i = 0; i < 2; ++i)
        for (unsigned int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(mtt(i,j), m(i,j));
}
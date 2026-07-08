// Core/DenseMatrix.h
#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

template <typename T>
class DenseMatrix {
private:
    unsigned int rows_;
    unsigned int cols_;
    T* data_;

    T Determinant2x2() const;
    T Determinant3x3() const;
    T DeterminantLU() const;
    DenseMatrix<T> Inverse1x1() const;
    DenseMatrix<T> Inverse2x2() const;
    DenseMatrix<T> Inverse3x3() const;
    DenseMatrix<T> InverseGaussJordan() const;

public:
    DenseMatrix() : rows_(0), cols_(0), data_(nullptr) {}
    DenseMatrix(unsigned int rows, unsigned int cols);
    DenseMatrix(const DenseMatrix<T>& mat);
    ~DenseMatrix() { delete[] data_; }

    DenseMatrix<T>& operator=(const DenseMatrix<T>& mat);
    void Resize(unsigned int rows, unsigned int cols);

    unsigned int GetRow () const { return rows_; }
    unsigned int GetCol () const { return cols_; }
    unsigned int GetSize() const { return rows_ * cols_; }

    T&       operator()(unsigned int row, unsigned int col);
    const T& operator()(unsigned int row, unsigned int col) const;

    void SetZero() { std::fill_n(data_, rows_ * cols_, T(0)); }

    T*       RawPtr()       { return data_; }
    const T* RawPtr() const { return data_; }

    std::vector<T>  DotVec  (const std::vector<T>& vec) const;
    DenseMatrix<T>  DotMat  (const DenseMatrix<T>& mat) const;
    DenseMatrix<T>  Transpose() const;
    T               Determinant() const;
    DenseMatrix<T>  Inverse() const;
};

// ======================= 实现 =======================

template <typename T>
DenseMatrix<T>::DenseMatrix(unsigned int rows, unsigned int cols)
    : rows_(rows), cols_(cols), data_(new T[rows * cols]()) {}

template <typename T>
DenseMatrix<T>::DenseMatrix(const DenseMatrix<T>& mat)
    : rows_(mat.rows_), cols_(mat.cols_), data_(new T[rows_ * cols_]) {
    std::copy(mat.data_, mat.data_ + rows_ * cols_, data_);
}

template <typename T>
DenseMatrix<T>& DenseMatrix<T>::operator=(const DenseMatrix<T>& mat) {
    if (this != &mat) {
        if (rows_ * cols_ != mat.rows_ * mat.cols_) {
            delete[] data_;
            data_ = new T[mat.rows_ * mat.cols_];
        }
        rows_ = mat.rows_;
        cols_ = mat.cols_;
        std::copy(mat.data_, mat.data_ + rows_ * cols_, data_);
    }
    return *this;
}

template <typename T>
void DenseMatrix<T>::Resize(unsigned int rows, unsigned int cols) {
    if (rows_ != rows || cols_ != cols) {
        delete[] data_;
        rows_ = rows;
        cols_ = cols;
        data_ = new T[rows * cols]();
    } else {
        SetZero();
    }
}

template <typename T>
T& DenseMatrix<T>::operator()(unsigned int row, unsigned int col) {
    assert(row < rows_ && col < cols_);
    return data_[col * rows_ + row];   // 列优先
}

template <typename T>
const T& DenseMatrix<T>::operator()(unsigned int row, unsigned int col) const {
    assert(row < rows_ && col < cols_);
    return data_[col * rows_ + row];
}

template <typename T>
std::vector<T> DenseMatrix<T>::DotVec(const std::vector<T>& vec) const {
    if (cols_ != vec.size())
        throw std::invalid_argument("Matrix-vector dimension mismatch");
    std::vector<T> out(rows_, T(0));
    for (unsigned int r = 0; r < rows_; ++r)
        for (unsigned int c = 0; c < cols_; ++c)
            out[r] += (*this)(r, c) * vec[c];
    return out;
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::DotMat(const DenseMatrix<T>& mat) const {
    if (cols_ != mat.rows_)
        throw std::invalid_argument("Matrix dimension mismatch");
    DenseMatrix<T> out(rows_, mat.cols_);
    out.SetZero();
    for (unsigned int j = 0; j < mat.cols_; ++j)
        for (unsigned int k = 0; k < cols_; ++k)
            for (unsigned int i = 0; i < rows_; ++i)
                out(i, j) += (*this)(i, k) * mat(k, j);
    return out;
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::Transpose() const {
    DenseMatrix<T> out(cols_, rows_);
    for (unsigned int r = 0; r < rows_; ++r)
        for (unsigned int c = 0; c < cols_; ++c)
            out(c, r) = (*this)(r, c);
    return out;
}

template <typename T>
T DenseMatrix<T>::Determinant() const {
    assert(rows_ == cols_);
    if      (rows_ == 1) return data_[0];
    else if (rows_ == 2) return Determinant2x2();
    else if (rows_ == 3) return Determinant3x3();
    else                 return DeterminantLU();
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::Inverse() const {
    assert(rows_ == cols_);
    if      (rows_ == 1) return Inverse1x1();
    else if (rows_ == 2) return Inverse2x2();
    else if (rows_ == 3) return Inverse3x3();
    else                 return InverseGaussJordan();
}

// ------- 私有方法实现 -------
template <typename T>
T DenseMatrix<T>::Determinant2x2() const {
    return (*this)(0,0) * (*this)(1,1) - (*this)(0,1) * (*this)(1,0);
}

template <typename T>
T DenseMatrix<T>::Determinant3x3() const {
    return (*this)(0,0) * ((*this)(1,1)*(*this)(2,2) - (*this)(1,2)*(*this)(2,1))
         - (*this)(0,1) * ((*this)(1,0)*(*this)(2,2) - (*this)(1,2)*(*this)(2,0))
         + (*this)(0,2) * ((*this)(1,0)*(*this)(2,1) - (*this)(1,1)*(*this)(2,0));
}

template <typename T>
T DenseMatrix<T>::DeterminantLU() const {
    DenseMatrix<T> A(*this);
    T det = T(1);
    for (unsigned int k = 0; k < rows_; ++k) {
        unsigned int pivot = k;
        T maxv = std::abs(A(k, k));
        for (unsigned int i = k + 1; i < rows_; ++i)
            if (std::abs(A(i, k)) > maxv) { maxv = std::abs(A(i, k)); pivot = i; }
        if (std::abs(A(pivot, k)) < T(1e-12)) return T(0);
        if (pivot != k) {
            for (unsigned int j = 0; j < cols_; ++j) std::swap(A(k, j), A(pivot, j));
            det = -det;
        }
        det *= A(k, k);
        for (unsigned int i = k + 1; i < rows_; ++i) {
            T f = A(i, k) / A(k, k);
            for (unsigned int j = k + 1; j < cols_; ++j)
                A(i, j) -= f * A(k, j);
        }
    }
    return det;
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::Inverse1x1() const {
    if (std::abs(data_[0]) < T(1e-12))
        throw std::runtime_error("Matrix is singular");
    DenseMatrix<T> r(1, 1);
    r(0, 0) = T(1) / data_[0];
    return r;
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::Inverse2x2() const {
    T det = Determinant2x2();
    if (std::abs(det) < T(1e-12))
        throw std::runtime_error("Matrix is singular");
    DenseMatrix<T> r(2, 2);
    r(0, 0) =  (*this)(1,1) / det;
    r(0, 1) = -(*this)(0,1) / det;
    r(1, 0) = -(*this)(1,0) / det;
    r(1, 1) =  (*this)(0,0) / det;
    return r;
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::Inverse3x3() const {
    T det = Determinant3x3();
    if (std::abs(det) < T(1e-12))
        throw std::runtime_error("Matrix is singular");
    DenseMatrix<T> r(3, 3);
    r(0,0)=((*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1))/det;
    r(0,1)=((*this)(0,2)*(*this)(2,1)-(*this)(0,1)*(*this)(2,2))/det;
    r(0,2)=((*this)(0,1)*(*this)(1,2)-(*this)(0,2)*(*this)(1,1))/det;
    r(1,0)=((*this)(1,2)*(*this)(2,0)-(*this)(1,0)*(*this)(2,2))/det;
    r(1,1)=((*this)(0,0)*(*this)(2,2)-(*this)(0,2)*(*this)(2,0))/det;
    r(1,2)=((*this)(0,2)*(*this)(1,0)-(*this)(0,0)*(*this)(1,2))/det;
    r(2,0)=((*this)(1,0)*(*this)(2,1)-(*this)(1,1)*(*this)(2,0))/det;
    r(2,1)=((*this)(0,1)*(*this)(2,0)-(*this)(0,0)*(*this)(2,1))/det;
    r(2,2)=((*this)(0,0)*(*this)(1,1)-(*this)(0,1)*(*this)(1,0))/det;
    return r;
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::InverseGaussJordan() const {
    DenseMatrix<T> aug(rows_, 2 * cols_);
    for (unsigned int i = 0; i < rows_; ++i) {
        for (unsigned int j = 0; j < cols_; ++j) {
            aug(i, j) = (*this)(i, j);
            aug(i, j + cols_) = (i == j) ? T(1) : T(0);
        }
    }
    for (unsigned int k = 0; k < rows_; ++k) {
        unsigned int pivot = k;
        T maxv = std::abs(aug(k, k));
        for (unsigned int i = k + 1; i < rows_; ++i)
            if (std::abs(aug(i, k)) > maxv) { maxv = std::abs(aug(i, k)); pivot = i; }
        if (std::abs(aug(pivot, k)) < T(1e-12))
            throw std::runtime_error("Matrix is singular");
        if (pivot != k)
            for (unsigned int j = 0; j < 2 * cols_; ++j)
                std::swap(aug(k, j), aug(pivot, j));
        T pv = aug(k, k);
        for (unsigned int j = 0; j < 2 * cols_; ++j) aug(k, j) /= pv;
        for (unsigned int i = 0; i < rows_; ++i) {
            if (i == k) continue;
            T f = aug(i, k);
            for (unsigned int j = 0; j < 2 * cols_; ++j)
                aug(i, j) -= f * aug(k, j);
        }
    }
    DenseMatrix<T> r(rows_, cols_);
    for (unsigned int i = 0; i < rows_; ++i)
        for (unsigned int j = 0; j < cols_; ++j)
            r(i, j) = aug(i, j + cols_);
    return r;
}
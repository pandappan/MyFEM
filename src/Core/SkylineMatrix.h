/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include <climits>
#include <vector>


template<class T_>
class DenseMatrix;

//! CSkylineMatrix class is used to store the FEM stiffness matrix in skyline storage
template <class T_>
class CSkylineMatrix
{
    
private:
//! Store the stiffness matrkix in skyline storage
    std::vector<T_> data_;
    
//! Dimension of the stiffness matrix
    unsigned int NEQ_;

//! Maximum half bandwith
    unsigned int MK_;

//! Size of the storage used to store the stiffness matrkix in skyline
    unsigned int NWK_;

//! Column hights
    std::vector<unsigned int> ColumnHeights_;
    
//! Diagonal address of all columns in data_
    std::vector<unsigned int> DiagonalAddress_;
    
public:

//! constructors
    inline CSkylineMatrix();
    inline CSkylineMatrix(unsigned int N);
    
//! destructor
    ~CSkylineMatrix() = default;

//! operator (i,j) where i and j numbering from 1
//! For the sake of efficiency, the index bounds are not checked
    inline T_& operator()(unsigned int i, unsigned int j);
    inline const T_& operator()(unsigned int i, unsigned int j) const;
    
//! Allocate storage for the skyline matrix
    inline void Allocate();
    
//! Calculate the column height, used with the skyline storage scheme
    void CalculateColumnHeight(const std::vector<unsigned int> &lm);

//! Calculate the maximum half bandwidth ( = max(ColumnHeights) + 1 )
    void CalculateMaximumHalfBandwidth();

//! Calculate address of diagonal elements in banded matrix
//! Caution: Address is numbered from 1 !
    void CalculateDiagnoalAddress();

//! Assemble the element stiffness matrix to the global stiffness matrix
    void Assembly(DenseMatrix<T_>& Matrix, const std::vector<unsigned int> &lm);

//! Return pointer to the ColumnHeights_
    inline std::vector<unsigned int>& GetColumnHeights();

//! Return the maximum half bandwidth
    inline unsigned int GetMaximumHalfBandwidth() const;

//! Return pointer to the DiagonalAddress_
    inline std::vector<unsigned int>& GetDiagonalAddress();

//! Return the dimension of the stiffness matrix
    inline unsigned int dim() const;
    
//! Return the size of the storage used to store the stiffness matrkix in skyline
    inline unsigned int size() const;

};

//! constructor functions
template <class T_>
inline CSkylineMatrix<T_>::CSkylineMatrix()
{
    NEQ_ = 0;
    MK_  = 0;
    NWK_ = 0;
}

template <class T_>
inline CSkylineMatrix<T_>::CSkylineMatrix(unsigned int N)
{
    NEQ_ = N;
    MK_  = 0;
    NWK_ = 0;

    ColumnHeights_.reserve(NEQ_);
    for (unsigned int i = 0; i < NEQ_; i++)
        ColumnHeights_[i] = 0;

    DiagonalAddress_.reserve(NEQ_ + 1);
    for (unsigned int i = 0; i < NEQ_ + 1; i++)
        DiagonalAddress_[i] = 0;
}

//! operator function (i,j) where i and j numbering from 1
template <class T_>
inline T_& CSkylineMatrix<T_>::operator()(unsigned int i, unsigned int j)
{
    if (j >= i)
        return data_[DiagonalAddress_[j - 1] + (j - i) - 1];
    else
        return data_[DiagonalAddress_[i - 1] + (i - j) - 1];
}

template <class T_>
inline const T_& CSkylineMatrix<T_>::operator()(unsigned int i, unsigned int j) const {
    if (j >= i)
        return data_[DiagonalAddress_[j - 1] + (j - i) - 1];
    else
        return data_[DiagonalAddress_[i - 1] + (i - j) - 1];
}

//! Allocate storage for the matrix
template <class T_>
inline void CSkylineMatrix<T_>::Allocate()
{
    NWK_ = DiagonalAddress_[NEQ_] - DiagonalAddress_[0];

    data_.assign(NWK_, 0.0);
    for (unsigned int i = 0; i < NWK_; i++)
        data_[i] = T_(0);
}

// 列高数组的引用
template <class T_>
inline std::vector<unsigned int>& CSkylineMatrix<T_>::GetColumnHeights()
{
    return ColumnHeights_;
}

//! Return the maximum half bandwidth
template <class T_>
inline unsigned int CSkylineMatrix<T_>::GetMaximumHalfBandwidth() const
{
    return(MK_);
}

//! Return pointer to the DiagonalAddress_
template <class T_>
inline std::vector<unsigned int>& CSkylineMatrix<T_>::GetDiagonalAddress()
{
    return DiagonalAddress_;
}

//! Return the dimension of the stiffness matrix
template <class T_>
inline unsigned int CSkylineMatrix<T_>::dim() const
{
    return(NEQ_);
}

//! Return the size of the storage used to store the stiffness matrkix in skyline
template <class T_>
inline unsigned int CSkylineMatrix<T_>::size() const
{
   return(NWK_);
}

//  Calculate the column height, used with the skyline storage scheme
template <class T_>
void CSkylineMatrix<T_>::CalculateColumnHeight(const std::vector<unsigned int>& lm)
{
    unsigned int nfirstrow = INT_MAX;
    for (unsigned int v : lm)
        if (v && v < nfirstrow) nfirstrow = v;
    for (unsigned int v : lm) {
        if (!v) continue;
        unsigned int h = v - nfirstrow;
        if (ColumnHeights_[v - 1] < h) ColumnHeights_[v - 1] = h;
    }
}

// Maximum half bandwidth ( = max(ColumnHeights) + 1 )
template <class T_>
void CSkylineMatrix<T_>::CalculateMaximumHalfBandwidth()
{
    MK_ = ColumnHeights_[0];

    for (unsigned int i=1; i<NEQ_; i++)
        if (MK_ < ColumnHeights_[i])
            MK_ = ColumnHeights_[i];

    MK_ = MK_ + 1;
}

//    Assemble the banded global stiffness matrix (skyline storage scheme)
template <class T_>
void CSkylineMatrix<T_>::Assembly(DenseMatrix<T_>& Matrix, const std::vector<unsigned int> &lm)
{
//  Assemble global stiffness matrix
    for (unsigned int j = 0; j < Matrix.GetCol(); j++) // 列循环
    {
        unsigned int Lj = lm[j];    // Global equation number corresponding to jth DOF of the element
        if (!Lj) continue; // 跳过约束自由度

        for (unsigned int i = 0; i <= j; i++) // 上三角部分
        {
            unsigned int Li = lm[i];    // Global equation number corresponding to ith DOF of the element

            if (!Li) continue; // 跳过约束自由度

            (*this)(Li,Lj) += Matrix(i, j);
        }
    }
}

//    Calculate address of diagonal elements in banded matrix
//    Caution: Address is numbered from 1 !
template <class T_>
void CSkylineMatrix<T_>::CalculateDiagnoalAddress()
{
    //    Calculate the address of diagonal elements
    //    M(0) = 1;  M(i+1) = M(i) + H(i) + 1 (i = 0:NEQ)
    DiagonalAddress_[0] = 1;
    for (unsigned int col = 1; col <= NEQ_; col++)
        DiagonalAddress_[col] = DiagonalAddress_[col - 1] + ColumnHeights_[col-1] + 1;

}


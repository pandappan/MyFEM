//
// Created by Administrator on 2026/6/11.
//
#include <iostream>
#include "Element.h"
#include "../Node.h"
#include "../../Core/Types.h"
#include "../../Core/DenseMatrix.h"

// 获取单元所有节点的坐标
// x_iI
CElement::CElement():ElementNumber_(0), NDim_(0), NEN_(0), ND_(0), ElementMaterial_(nullptr), volume_((0.0)){}

DenseMatrix<double> CElement::GetNodeCoordinates() const {
    DenseMatrix<double> nodeCoords(NDim_, NEN_);
    for (unsigned int k = 0; k < NEN_; k++) {
        for (unsigned int d = 0; d < NDim_; d++) {
            nodeCoords(d, k) = nodes_[k]->XYZ[d];
        }
    }
    return nodeCoords;
}

void CElement::GenerateLocationMatrix()
{
    const DOFIndex* activeDOFs  = GetActiveDOFs();
    unsigned int ndofs = GetNumActiveDOFsPerNode();
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < ndofs; D++)
            LocationMatrix_[i++] = nodes_[N]->eqn[activeDOFs[D]];
}

void CElement::AllocateStorage(unsigned int nDim, unsigned int nen, unsigned int nd) {
    NEN_ = nen;
    NDim_ = nDim;
    ND_ = nd;
    nodes_.assign(NEN_, nullptr);
    LocationMatrix_.assign(ND_, 0);
}
void CElement::SetUpForTesting(std::vector<CNode*>& NodeList, CMaterial* Material_) {
    nodes_ = std::move(NodeList);
    ElementMaterial_ = Material_;
}
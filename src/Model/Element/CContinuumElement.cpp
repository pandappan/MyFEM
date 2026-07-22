//
// Created by Administrator on 2026/6/10.
//

#include <cmath>
#include <stdexcept>
#include "../Node.h"
#include "CContinuumElement.h"
#include "../Material/Material.h"

// 计算雅可比矩阵
// Jac_ij = (ax/axi) = x_iI o (aN_I/axi)_jI
void CContinuumElement::ComputeJacobian(const DenseMatrix<double>& dN_dxi,
    const DenseMatrix<double> &nodeCoords, DenseMatrix<double> &jacobian) const {
    jacobian.Resize(NDim_, NDim_);
    jacobian.SetZero();
    for (unsigned int i = 0; i < NDim_; i++)
        for (unsigned int j = 0; j < NDim_; j++)
            for (unsigned int I = 0; I < NEN_; I++)
                jacobian(i, j) += nodeCoords(i, I) * dN_dxi(j, I);
}

// 计算雅可比矩阵的逆以及行列式
// (Jac^-1)_ij detJac
double CContinuumElement::ComputeInverseJacobian(const DenseMatrix<double> &jacobian,
    DenseMatrix<double> &invJacobian) const {
    double detJ = jacobian.Determinant();
    if (std::abs(detJ) < 1e-12) {
        throw std::runtime_error("Jacobian determinant is zero or near zero");
    }
    invJacobian = jacobian.Inverse();
    return detJ;
}

// 计算形函数的全局导数
// dN_I/dx_i = dN_I/dxi_j * (Jac^-1)_ji
void CContinuumElement::ComputeGlobalDerivatives(const DenseMatrix<double>& dN_dxi,
                              const DenseMatrix<double>& invJacobian,
                              DenseMatrix<double>& dN_dx) const {
    dN_dx.SetZero();
    for (unsigned int i = 0; i < NDim_; i++) {
        for (unsigned int j = 0; j < NDim_; j++) {
            for (unsigned int k = 0; k < NEN_; k++) {
                dN_dx(i, k) +=  dN_dxi(j,k) * invJacobian(j,i);
            }
        }
    }
}

// 预先计算一个积分点处所有信息，在所有积分点处循环调用，进而形成积分点处的所有信息
void CContinuumElement::ComputeIntegrationPointData(
    const std::vector<double>& xi,
    double weight,
    const DenseMatrix<double>& nodeCoords,
    IntegrationPointData& ipData) {
    // 计算形函数值
    ComputeShapeFunctions(xi, ipData.N);
    // 计算形函数局部导数
    DenseMatrix<double> dN_dxi(NDim_, NEN_);
    ComputeShapeDerivatives(xi, dN_dxi);
    // 计算雅可比矩阵
    DenseMatrix<double> jacobian(NDim_, NDim_);
    ComputeJacobian(dN_dxi, nodeCoords, jacobian);
    // 计算雅可比逆矩阵和行列式
    DenseMatrix<double> invJacobian(NDim_, NDim_);
    double detJ = ComputeInverseJacobian(jacobian, invJacobian);
    // 计算形函数全局导数
    ComputeGlobalDerivatives(dN_dxi, invJacobian, ipData.dN_dx);
    // 计算 detJ * weight
    ipData.detJ_times_weight = detJ * weight;
}

// 计算B矩阵
void CContinuumElement::ComputeBMatrix(unsigned int ip,
    DenseMatrix<double>& B) const {
    // 引用语法
    const auto& dN_dx = integrationPoints_[ip].dN_dx;
    if (NDim_ == 2) {
        for (unsigned int i = 0; i < NEN_; i++) {
            double dNdx = dN_dx(0,i);
            double dNdy = dN_dx(1, i);
            B(0, 2*i    ) = dNdx;
            B(1, 2*i + 1) = dNdy;
            B(2, 2*i    ) = dNdy;
            B(2, 2*i + 1) = dNdx;
        }
    } else if (NDim_ == 3) {
        for (unsigned int i = 0; i < NEN_; i++) {
            double dNdx = dN_dx(0, i);
            double dNdy = dN_dx(1, i);
            double dNdz = dN_dx(2, i);
            B(0, 3*i    ) = dNdx;
            B(1, 3*i + 1) = dNdy;
            B(2, 3*i + 2) = dNdz;
            B(3, 3*i    ) = dNdy;  B(3, 3*i + 1) = dNdx;
            B(4, 3*i + 1) = dNdz;  B(4, 3*i + 2) = dNdy;
            B(5, 3*i    ) = dNdz;  B(5, 3*i + 2) = dNdx;
        }
    }
}

// 循环调用ComputeIntegrationPointData形成单元所有积分点处的信息，
// 包括形函数，全局导数，雅可比行列*权重，体积
// 初始一次性形成，后续计算刚度矩阵，体力，约束力时直接使用
void CContinuumElement::InitializeIntegrationPoints() {
    if (integrationPointsCached_) {
        return;
    }
    // 积分点个数
    unsigned int nGp = GetNumIntegrationPoints();
    // 获取积分点位置以及权重
    GaussData Gauss = GetIntegrationPoint();
    // 调整积分点容器大小
    integrationPoints_.resize(nGp);
    // 为每个积分点初始化数据结构
    for (unsigned int i = 0; i < nGp; i++) {
        integrationPoints_[i].N.resize(NEN_);
        integrationPoints_[i].dN_dx.Resize(NDim_,NEN_);
        integrationPoints_[i].detJ_times_weight = 0.0;
    }
    // 单元的节点坐标数组
    DenseMatrix<double> nodeCoords = GetNodeCoordinates();
    // 计算积分点处的形函数信息
    for (unsigned int i = 0; i < nGp; i++) {
        // 提取积分点处的坐标
        std::vector<double> xi(NDim_);
        for (unsigned int j = 0; j < NDim_; j++) {
            xi[j] = Gauss.GaussPoints(i, j);
        }
        // 计算该积分点的所有信息
        ComputeIntegrationPointData(xi, Gauss.GaussWeights[i],
                                    nodeCoords, integrationPoints_[i]);
    }
    volume_ = 0.0;
    CMaterial* mat = GetElementMaterial();
    double thk = 1.0;
    if (mat) { // 单元正常作为体元的分支
        for (unsigned int ip = 0; ip < nGp; ip++) {
            thk = GetIntegrationVolumeFactor(ip);
            volume_ += integrationPoints_[ip].detJ_times_weight * thk;
        }
    } else { // 单元作为载荷面元的分支
        for (unsigned int ip = 0; ip < nGp; ip++) {
            volume_ += integrationPoints_[ip].detJ_times_weight * thk;
        }
    }
    integrationPointsCached_ = true;
}

// Ke = B^T o D o B * detJac * W_I
void CContinuumElement::ElementStiffness(DenseMatrix<double> &Ke) const {
    Ke.SetZero();
    CMaterial* mat = GetElementMaterial();
    unsigned int numStrain = GetNumStrainComponents();
    DenseMatrix<double> D(numStrain, numStrain);
    mat->ComputeElasticMatrix(D);
    DenseMatrix<double> B(numStrain, ND_);
    for (unsigned int ip = 0; ip < GetNumIntegrationPoints(); ++ip) {
        ComputeBMatrix(ip, B);
        double thk = GetIntegrationVolumeFactor(ip);
        double dv = integrationPoints_[ip].detJ_times_weight * thk;
        DenseMatrix<double> DB   = D.DotMat(B);              // (nS, ND_)
        DenseMatrix<double> BtDB = B.Transpose().DotMat(DB); // (ND_, ND_)
        for (unsigned int i = 0; i < ND_; ++i)
            for (unsigned int j = 0; j < ND_; ++j)
                Ke(i, j) += BtDB(i, j) * dv;
    }
}

// 体载转化为点载荷并写入节点中
// Fext_I = N_I o b * rho * detJac * W_J
void CContinuumElement::CalculateBodyForce(const double *bodyForce) {
    // 单元的基本材料参数
    CMaterial* mat = GetElementMaterial();
    double rho = mat->rho;
    if (std::abs(rho) < 1.0e-12) return;
    // 积分点缓存变量
    if (!integrationPointsCached_) InitializeIntegrationPoints();
    for (unsigned int i = 0; i < NEN_; i++) {
        // 累加计算系数
        double factor = 0.0;
        for (unsigned int ip = 0; ip < GetNumIntegrationPoints(); ip++) {
            double thk = GetIntegrationVolumeFactor(ip);
            double dv = integrationPoints_[ip].detJ_times_weight * thk;
            factor += integrationPoints_[ip].N[i] * rho * dv;
        }
        // 计算等效节点力，并存入节点
        for (unsigned int d = 0; d < NDim_; d++) {
            double eqforce = factor * bodyForce[d];
            nodes_[i]->AddForce(d,eqforce);
        }
    }
}

// epsilon = B o u
std::vector<double> CContinuumElement::ComputeStrainAtIntegrationPoint(unsigned int ip) const {
    unsigned int numStrain = GetNumStrainComponents();
    DenseMatrix<double> B(numStrain, ND_);
    ComputeBMatrix(ip, B);
    std::vector<double> ue(ND_, 0.0);
    const DOFIndex* dofs = GetActiveDOFs();
    unsigned int ndofs = GetNumActiveDOFsPerNode();
    int index = 0;
    for (unsigned int i = 0; i < NEN_; i++) {
        for (unsigned int j = 0; j < ndofs; j++) {
            ue[index++] = nodes_[i]->displacement[dofs[j]];
        }
    }
    std::vector<double> strain(numStrain, 0.0);
    strain = B.DotVec(ue);
    return strain;
}

std::vector<double> CContinuumElement::ComputeStressAtIntegrationPoint(unsigned int ip) const {
    unsigned int numStrain = GetNumStrainComponents();
    std::vector<double> strain(numStrain, 0.0);
    strain = ComputeStrainAtIntegrationPoint(ip);
    std::vector<double> stress(numStrain, 0.0);
    GetElementMaterial()->ComputeStress(strain, stress);
    return stress;
}

// (NEN_ x nGp)，默认简单平均
DenseMatrix<double> CContinuumElement::GetExprapolationMatrix() const {
    const unsigned int nGp = GetNumIntegrationPoints();
    DenseMatrix<double> E(NEN_, nGp);
    for (unsigned int i = 0; i < NEN_; i++) {
        for (unsigned int j = 0; j < nGp; j++) {
            E(i,j) = 1.0 / nGp;
        }
    }
    return E;
}

void CContinuumElement::ExtrapolatStressToNodes(std::vector<std::vector<double> > &nodalStress) const {
    // 积分点个数，应力分量数目
    const unsigned int nGp = GetNumIntegrationPoints();
    const unsigned int numStrain = GetNumStrainComponents();
    // 所有积分点处的应力
    std::vector<std::vector<double>> gpStress(nGp);
    for (unsigned int i = 0; i < nGp; i++) {
        gpStress[i] = ComputeStressAtIntegrationPoint(i);
    }
    // 外推矩阵
    DenseMatrix<double> E = GetExprapolationMatrix();
    // 初始化输出
    nodalStress.assign(NEN_, std::vector<double>(numStrain, 0.0));
    // 逐节点，逐分量外推应力
    for (unsigned int n = 0; n < NEN_; n++) {
        for (unsigned int c = 0; c < numStrain; c++) {
            double stress = 0.0;
            for (unsigned int gp = 0; gp < nGp; gp++) {
                stress += E(n,gp) * gpStress[gp][c];
            }
            nodalStress[n][c] = stress;
        }
    }
}

// 返回单元所有积分点处的位置
DenseMatrix<double> CContinuumElement::GetIntegrationPointPositions() const {
    const unsigned int nGp = GetNumIntegrationPoints();
    DenseMatrix<double> pos(NDim_, nGp);
    pos.SetZero();
    for (unsigned int ip = 0; ip < nGp; ip++) {
        const auto& N = integrationPoints_[ip].N;
        for (unsigned int d = 0; d < NDim_; d++) {
            for (unsigned int n = 0; n < NEN_; n++) {
                pos(d,ip) += N[n] * nodes_[n]->XYZ[d];
            }
        }
    }
    return pos;
}

// 返回所有积分点处的应力(nGp,nComp)
std::vector<std::vector<double> > CContinuumElement::GetIntegrationPointStresses() const {
    const unsigned int nGp = GetNumIntegrationPoints();
    const unsigned int numStrain = GetNumStrainComponents();
    std::vector<std::vector<double>> stress(nGp, std::vector<double>(numStrain, 0.0));
    for (unsigned int i = 0; i < nGp; i++) {
        stress[i] = ComputeStressAtIntegrationPoint(i);
    }
    return stress;
}

double CContinuumElement::GetRepresentativeStress() const {
    // 将所有积分点的miss应力平均
    unsigned int nGp = GetNumIntegrationPoints();
    unsigned int numStrain = GetNumStrainComponents();
    std::vector<std::vector<double>> stress = GetIntegrationPointStresses();
    double miss = 0.0;
    for (const auto& ipStress: stress) {
        double temp = 0.0;
        if (numStrain == 3) {
            double sxx = ipStress[0];
            double syy = ipStress[1];
            double sxy = ipStress[2];
            temp = std::sqrt(sxx*sxx + syy*syy - sxx*syy + 3.0 * sxy*sxy);

        } else if (numStrain == 6) {
            double sxx=ipStress[0], syy=ipStress[1], szz=ipStress[2];
            double sxy=ipStress[3], syz=ipStress[4], sxz=ipStress[5];
            double d1 = sxx-syy, d2 = syy-szz, d3 = szz-sxx;
            temp = std::sqrt(0.5*(d1*d1+d2*d2+d3*d3)
                             + 3.0*(sxy*sxy+syz*syz+sxz*sxz));
        }
        miss += temp;
    }
    return miss / nGp;
}

// r = sum N_I xi_I
double CContinuumElement::GetRadiusAtIntegrationPoint(unsigned int ip) const {
    const std::vector<double> N = integrationPoints_[ip].N;
    double radius = 0.0;
    for (unsigned int I = 0; I < NEN_; I++) {
        radius += N[I] * nodes_[I]->XYZ[0];
    }
    return radius;
}

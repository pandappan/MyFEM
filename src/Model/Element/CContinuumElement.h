//
// Created by Administrator on 2026/6/10.
//

#pragma once

#include <vector>
#include "Element.h"
#include "../../Core/DenseMatrix.h"

template<class T>
class DenseMatrix;

// 积分点信息
struct GaussData {
    std::vector<double> GaussWeights;
    DenseMatrix<double> GaussPoints;
};
// 连续介质单元基类
class CContinuumElement : public CElement {
protected:
    struct IntegrationPointData {
        // 形函数
        std::vector<double> N;
        // 形函数的全局导数
        DenseMatrix<double> dN_dx;
        // 雅可比矩阵*权重
        double detJ_times_weight;
    };
    std::vector<IntegrationPointData> integrationPoints_;
    bool integrationPointsCached_ = false;
    // 单元的外推矩阵，如果没有特殊外推矩阵，则取各节点取积分平均值
    virtual DenseMatrix<double> GetExprapolationMatrix() const;
public:
    //========通用逻辑========
    // 单元刚度矩阵
    void ElementStiffness(DenseMatrix<double>& Ke) const override;
    // 连续介质单元独有逻辑
    void OnSetupComplete() override {InitializeIntegrationPoints();};
    // 计算所有积分点的形函数信息
    void InitializeIntegrationPoints();
    // 将体载转化为点载
    void CalculateBodyForce(const double *bodyForce) override;
    // 计算积分点处的应变
    std::vector<double> ComputeStrainAtIntegrationPoint(unsigned int ip) const;
    // 计算单元在积分点处的应力
    std::vector<double> ComputeStressAtIntegrationPoint(unsigned int ip) const;
    // 单元平均的miss应力
    double GetRepresentativeStress() const override;

    //========查询函数========
    // 获取单元积分点的数目
    const std::vector<IntegrationPointData>& GetIntegrationPoints() const {
        return integrationPoints_;
    }

    //========派生逻辑========
    virtual unsigned int GetNumIntegrationPoints() const = 0;
    virtual GaussData GetIntegrationPoint() const = 0;
    virtual void ComputeShapeFunctions(const std::vector<double>& xi,
        std::vector<double>& N) const = 0;
    virtual void ComputeShapeDerivatives(const std::vector<double>& xi,
        DenseMatrix<double>& dN_dxi) const = 0;
    DenseMatrix<double> GetIntegrationPointPositions() const;
    std::vector<std::vector<double>> GetIntegrationPointStresses() const;
    // 将所有积分点的应力外推至节点应力
    virtual void ExtrapolatStressToNodes(std::vector<std::vector<double>>& nodalStress) const;

private:
    //========私有辅助========
    //! 计算应变B矩阵
    void ComputeBMatrix(unsigned int ip, DenseMatrix<double>& BMatrix) const;
    //! 计算雅可比矩阵
    void ComputeJacobian(const DenseMatrix<double>& dN_dxi,
        const DenseMatrix<double>& nodeCoords,
        DenseMatrix<double>& jacobian) const;
    //! 计算雅可比矩阵的逆
    double ComputeInverseJacobian(const DenseMatrix<double>& jacobian,
        DenseMatrix<double>& invJacobian) const;
    //! 计算形函数的全局导数
    void ComputeGlobalDerivatives(const DenseMatrix<double>& dN_dxi,
        const DenseMatrix<double>& invJacobian,
        DenseMatrix<double>& dN_dx) const;
    //! 计算积分点处的形函数等信息
    void ComputeIntegrationPointData(const std::vector<double>& xi,
        double weight,
        const DenseMatrix<double>& nodeCoords,
        IntegrationPointData& ipData);
};


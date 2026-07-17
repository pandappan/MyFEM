//
// Created by Administrator on 2026/6/10.
//
#pragma once
#include "CContinuumElement.h"
#include <vector>

class CQ4 : public CContinuumElement{
private:
    // 静态成员变量，属于类本身的特性，不单独属于某个类，静态常量可以在类内定义，静态常量数组需要在类外定义
    static const DOFIndex ActiveDOFs[2];
    static const unsigned int NumActiveDOFsPerNode;
public:
    // 构造函数
    CQ4();

    //! Q4单元作为H8单元的面元时接受节点数组初始化
    void AsFaceElem(const std::vector<CNode*>& nodelist);

    //!	Write element data to stream
    void Write(std::ostream& output) const override;
    void WriteElementStress(std::ostream& output) const override;
    //! 获取节点激活自由度数目
    unsigned int GetNumActiveDOFsPerNode() const override {return NumActiveDOFsPerNode;};

    //! 获取激活自由度局部编号
    const DOFIndex* GetActiveDOFs() const override {return ActiveDOFs;};

    //! 输出单元类型
    std::string ElementTypeName() const override {return "CQ4";};

    //! 返回节点初始坐标
    void GetVisualizationNodes(DenseMatrix<double>& coords) const;

    //! 返回节点变形坐标
    void GetVisualizationDeformeNodes(DenseMatrix<double>& deformeCoords) const;

    //! 积分点个数
    unsigned int GetNumIntegrationPoints() const override {
        return 4;
    }
    //! 高斯点的数据信息
    GaussData GetIntegrationPoint() const override;
    //! 形函数
    void ComputeShapeFunctions(const std::vector<double>& xi,
                          std::vector<double>& N) const override;
    //! 形函数的局部导数
    void ComputeShapeDerivatives(const std::vector<double>& xi,
                            DenseMatrix<double>& dN_dxi) const override;
    //! 面载转化为点载
    bool CalculateSurfaceLoad(unsigned int faceID, unsigned int dof, double value) override;
    //! 获取局部节点编号
    std::vector<int> GetFaceNodesLocalID(unsigned int faceID) const;
    //! 单元外推矩阵
    DenseMatrix<double> GetExprapolationMatrix() const override;
};

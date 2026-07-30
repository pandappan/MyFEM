#pragma once

#include <fstream>
#include <ostream>
#include <string>
#include <vector>
#include "../../Core/Types.h"

class Node;
class Material;
class Writer;
class CElementGroup;
template<class T>
class DenseMatrix;

//!	Element base class
class Element
{
protected:
    ElementTypes elementType_;
    unsigned int ElementNumber_;
    unsigned int NDim_; // 单元坐标维度
    unsigned int NEN_; // 单元节点数目
    unsigned int ND_;
    std::vector<Node*> nodes_;
    Material* ElementMaterial_;
    std::vector<unsigned int> LocationMatrix_;
    double volume_; // 单元体积

public:

    Element();
    virtual ~Element() = default;
    // 纯粹虚接口
    // 设置单元基础信息
    void SetElementInfo(unsigned int elemId_0, Material* matPtr,
        const std::vector<unsigned int>& connectivity_0, std::vector<Node>& nodeList);
    // 写出单元数据
    virtual void Write(std::ostream& output) const = 0;
    virtual void WriteElementStress(std::ostream& output) const = 0;
    // 计算单元刚度矩阵
    virtual void ElementStiffness(DenseMatrix<double>& Ke) const = 0;
    // 节点自由度数目
    virtual unsigned int GetNumActiveDOFsPerNode() const = 0;
    // 节点自由度编号
    virtual const DOFIndex* GetActiveDOFs() const = 0;
    // 单元类型
    virtual std::string ElementTypeName() const = 0;
    // 可视化接口
    // 返回可视化的节点列表
    virtual void GetVisualizationNodes(DenseMatrix<double>& coords) const = 0;
    // 返回节点位移场
    virtual void GetVisualizationDeformeNodes(DenseMatrix<double>& deformeCoords) const = 0;
    // 生成单元定位数组
    void GenerateLocationMatrix();
    // 单元编号
    void SetElementNumber(unsigned int n) {ElementNumber_ = n;}
    unsigned int GetElementNumber() const {return ElementNumber_;}
    // 单元节点数目
    unsigned int GetNEN() const {return NEN_;}
    // 单元节点指针数组
    const std::vector<Node*>& GetNodes() const {return nodes_;};
    // 单元材料指针
    Material* GetElementMaterial() const {return ElementMaterial_;};
    // 单元定位数组
    const std::vector<unsigned int>& GetLocationMatrix() const {return LocationMatrix_;};
    // 单元自由度总数
    unsigned int GetND() const {return ND_;};
    // 单元体积
    double GetVolume() const {return volume_;};
    // 获取单元所有节点的坐标
    DenseMatrix<double> GetNodeCoordinates() const;
    // 获取单元从属节点的位移以及约束状态
    void GetElementNodesDisp(std::vector<double>& nodesDisp, std::vector<int>& nodesBcode) const;
    // 获取单元从属节点的外载荷
    void GetElementNodesForce(std::vector<double>& nodesForce);
    // 单元层面回代求解反力，写入节点
    void CalculateBCForce();
    // 将面载荷转化为节点载荷，并存入节点中
    virtual bool CalculateSurfaceLoad(unsigned int faceID, unsigned int dof, double value);
    // 设置单元的类型
    void SetElementType(ElementTypes elementType);
    // 将体载荷转化为点载荷，并存入节点中
    virtual void CalculateBodyForce(const double* bodyForce);
    // 单元平均应力
    virtual double GetRepresentativeStress() const {return 0.0;};
    // 计算单元总应变能
    double CalculateElementEnergy() const;
    // 测试使用
    void SetupForTesting(std::vector<Node*> NodeList, Material* Material_);
    // 连续介质单元积分点信息完成初始化
    virtual void OnSetupComplete() {};
    // 声明单元所需的材料类型，用于兼容性检查
    virtual MaterialCategory GetRequiredMaterial() const = 0;
    // 将单元使用的自由度在所有节点上的掩码置位
    void RegisterDofsOnNodes();
protected:
    void AllocateStorage(unsigned int nDim, unsigned int nen, unsigned int nd);
};
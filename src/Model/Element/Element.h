#pragma once

#include <fstream>
#include <ostream>
#include <string>
#include <vector>
#include "../../Core/Types.h"

class CNode;
class CMaterial;
class Writer;
class CElementGroup;
template<class T>
class DenseMatrix;

//!	Element base class
class CElement
{
protected:
    unsigned int ElementNumber_;
    unsigned int NDim_; // 单元坐标维度
    unsigned int NEN_; // 单元节点数目
    unsigned int ND_;
    std::vector<CNode*> nodes_;
    CMaterial* ElementMaterial_;
    std::vector<unsigned int> LocationMatrix_;
    double volume_; // 单元体积

public:

    CElement();
    virtual ~CElement() = default;
    // 纯粹虚接口
    // 读取单元数据
    virtual bool Read(std::ifstream& Input, CElementGroup& Group, std::vector<CNode>& NodeList) = 0;
    // 写出单元数据
    virtual void Write(std::ostream& output) const = 0;
    virtual void WriteElementStress(std::ostream& output) const = 0;
    // 计算单元刚度矩阵
    virtual void ElementStiffness(DenseMatrix<double>& K) = 0;
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
    const std::vector<CNode*>& GetNodes() const {return nodes_;};
    // 单元材料指针
    CMaterial* GetElementMaterial() const {return ElementMaterial_;};
    // 单元定位数组
    const std::vector<unsigned int>& GetLocationMatrix() const {return LocationMatrix_;};
    // 单元自由度总数
    unsigned int GetND() const {return ND_;};
    // 单元体积
    double GetVolume() const {return volume_;};
    // 获取单元所有节点的坐标
    DenseMatrix<double> GetNodeCoordinates() const;
    // 测试使用
    void SetupForTesting(std::vector<CNode*>& NodeList, CMaterial* Material_);
protected:
    void AllocateStorage(unsigned int nDim, unsigned int nen, unsigned int nd);
};
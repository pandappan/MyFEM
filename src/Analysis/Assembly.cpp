//
// Created by Administrator on 2026/7/6.
//
#include <iostream>
#include "Assembly.h"
#include "../Model/Element/Element.h"
#include "../Model/Model.h"
#include "../Core/DenseMatrix.h"
#include "Element/CContinuumElement.h"
#include "Material/Material.h"

void Assembler::CalculateEquationNumber(Model& model) {
    model.neq = 0;
    for (auto& node : model.nodes)
        node.GenerateNodeEquation(model.neq);
}

void Assembler::CalculateLocationMatrix(Model &model) {
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            element.GenerateLocationMatrix();
        }
    }
}

void Assembler::AllocateLinearSystem(Model &model) {
    // Allocate for global force/displacement vector
    model.force.assign(model.neq, 0.0);

    // Create the banded stiffness matrix
    model.K.reset(new CSkylineMatrix<double>(model.neq));
    // 计算列高
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            model.K->CalculateColumnHeight(element.GetLocationMatrix());
        }
    }
    model.K->CalculateMaximumHalfBandwidth();
    model.K->CalculateDiagnoalAddress();
    model.K->Allocate();
}

// 填充单元映射表值：单元编号-单元指针，为后续面元寻找对应单元提供索引
void Assembler::InitializeElementMap(Model &model) {
    // 初始化内存
    unsigned int totolElems = 0;
    for (auto& group : model.groups) {
        totolElems += group.GetNUME();
    }
    model.globalElementMap.clear();
    model.globalElementMap.reserve(totolElems);
    // 填充映射表值
    for (auto& group : model.groups) {
        for (unsigned int e = 0; e < group.GetNUME(); e++) {
            CElement* elem = &(group.GetElement(e));
            unsigned int elemId_0 = elem->GetElementNumber();
            auto it = model.globalElementMap.find(elemId_0);
            if (it != model.globalElementMap.end()) {
                throw std::runtime_error("Error: Element \"" + std::to_string(elemId_0) + "\" already exists");
            }
            model.globalElementMap.emplace(elemId_0, elem);
        }
    }
}

// 面载荷转化为等效节点力，存入节点中
void Assembler::ConvertSLoadsToCLoads(Model &model) {
    // 循环所有面载荷，逐步转化为等效节点载荷，并存入节点中
    for (auto& sload: model.sloads) {
        CElement* elem = model.globalElementMap[sload.elemId_0];
        elem->CalculateSurfaceLoad(sload.faceId_0, sload.dof_0, sload.value);
    }
}

// 所有单元体载荷转化为等效节点力，存入节点中
void Assembler::ConvertBLoadsToCLoads(Model &model) {
    // 体力加速度为0，直接退出，避免后续浪费计算量
    const double* b = model.bodyForce;
    double totalForce = 0.0;
    totalForce = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    if (totalForce < 1.0e-12) return;
    // 正常转换体力
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            element.CalculateBodyForce(b);
        }
    }
}

// 装配节点上的所有力，包括点载荷，面载的等效点载，体载的等效点载
void Assembler::AssembleForce(Model &model) {
    std::fill(model.force.begin(), model.force.end(), 0.0);
    for (auto& node : model.nodes) {
        for (unsigned int d = 0; d < CNode::NDF; d++) {
            unsigned int eq = node.eqn[d];
            if (eq) {
                model.force[eq - 1] += node.GetForce(d);
            }
        }
    }
}

// 装配单元刚度矩阵
// 计算并且装配位移约束的造成的右端修正项
// 受主从自由度影响，将单元刚度矩阵和右端项进行修正
void Assembler::AssembleStiffnessAndConstraintCorrection(Model &model) {
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            unsigned int nd = element.GetND();
            const std::vector<unsigned int>& lm = element.GetLocationMatrix();
            // 装配单元刚度矩阵
            DenseMatrix<double> ke(nd, nd);
            element.ElementStiffness(ke);
            model.K->Assembly(ke,lm);
            // 装配右端修正项
            std::vector<double> right(nd);
            element.ElementRight(ke, right);
            for (unsigned int i = 0; i < nd; i++) {
                if (lm[i] != 0) {
                    model.force[lm[i]-1] -= right[i];
                }
            }
        }
    }
}

void Assembler::WriteDisplacementToNodes(Model &model) {
    for (auto& node : model.nodes) {
        node.UpdataNodeDisplacement(model.force);
    }
}

void Assembler::CalculateNodalBCForce(Model &model) {
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            element.CalculateBCForce();
        }
    }
}

// 外推积分点应力，并进行面积加权平均。杆单元不需要外推到节点应力，只有连续介质单元需要处理
void Assembler::CalculateNodalStress(Model& model) {
    // 应力分量数目+判定单元类型
    unsigned int nComp = 0;
    for (auto& g : model.groups) {
        for (unsigned int e = 0; e < g.GetNUME(); e++) {
            if (g.GetNUME() > 0) {
                const auto* c = dynamic_cast<const CContinuumElement*>(&g.GetElement(e));
                if (c) {
                    nComp = c->GetElementMaterial()->GetNumStressComponents();
                    break;
                }
            }
        }
        if (nComp > 0) break;
    }
    // 桁架单元，跳过外推
    if (nComp == 0) return;
    // 初始化节点的应力容器
    for (auto& node : model.nodes) {
        node.stress.assign(nComp, 0.0);
        node.stressWieghts = 0.0;
    }
    // 遍历所有单元，外推+加权累加到节点
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            // 只有连续介质单元才需要外推
            auto* continuum = dynamic_cast<CContinuumElement*>(&element);
            if (!continuum) continue;
            // 面积权重
            double weight = continuum->GetVolume();
            if (weight <= 0.0) continue;
            // 外推
            std::vector<std::vector<double>> nodalStress;
            continuum->ExtrapolatStressToNodes(nodalStress);
            // 累加到全局节点
            const auto& elemNodes = continuum->GetNodes();
            for (unsigned int n = 0; n < elemNodes.size(); n++) {
                CNode* node = elemNodes[n];
                // 应力*权重
                for (unsigned int c = 0; c < nComp; c++) {
                    node->stress[c] += weight * nodalStress[n][c];
                }
                // 节点累计权重
                node->stressWieghts += weight;
            }
        }
    }
    // 累计加权应力/总权重得到节点平均应力
    for (auto& node: model.nodes) {
        for (unsigned int c = 0; c < nComp; c++) {
            node.stress[c] /= node.stressWieghts;
        }
    }
}
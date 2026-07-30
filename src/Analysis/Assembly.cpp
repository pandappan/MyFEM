//
// Created by Administrator on 2026/7/6.
//
#include <iostream>
#include <stdexcept>
#include "Assembly.h"
#include "../Model/Element/Element.h"
#include "../Model/Model.h"
#include "../Core/DenseMatrix.h"
#include "Element/ContinuumElement.h"
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
            Element& element = group.GetElement(e);
            element.GenerateLocationMatrix();
        }
    }
}

void Assembler::AllocateLinearSystem(Model &model) {
    // Allocate for global force/displacement vector
    model.force.assign(model.neq, 0.0);

    // Create the banded stiffness matrix
    model.K.reset(new SkylineMatrix<double>(model.neq));
    // 计算列高
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            Element& element = group.GetElement(e);
            std::vector<unsigned int> eqs = GetEffectiveEquations(element, model);
            model.K->CalculateColumnHeight(eqs);
        }
    }
    model.K->CalculateMaximumHalfBandwidth();
    model.K->Diagonal();
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
            Element* elem = &(group.GetElement(e));
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
        Element* elem = model.globalElementMap[sload.elemId_0];
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
            Element& element = group.GetElement(e);
            element.CalculateBodyForce(b);
        }
    }
}

// 装配节点外载荷到方程右端
// 自由度直接进方程；从自由度上的力按 T 分摊到 master
void Assembler::AssembleForce(Model &model) {
    std::fill(model.force.begin(), model.force.end(), 0.0);
    for (auto& node : model.nodes) {
        for (unsigned int d = 0; d < NDF_MAX; d++) {
            if (!node.IsDofActive(d)) continue; // 跳过未激活自由度
            double f = node.GetForce(d);
            if (f == 0.0) continue;

            switch (node.bcode[d]) {
                case 0:  // 自由：直接进方程
                    model.force[node.eqn[d] - 1] += f;
                    break;
                case 3: { // 从自由度：按约束关系分摊到 master
                    int m = model.FindMPCBySlave(node.Index, d);
                    if (m < 0) break;
                    for (const auto& t : model.mpcs[m].masters) {
                        const Node& mn = model.nodes[t.node_0];
                        if (mn.bcode[t.dof_0] == 0)   // 只有自由 master 进方程
                            model.force[mn.eqn[t.dof_0] - 1] += t.coeff * f;
                        // master 固定/指定位移：力落到约束上，算作反力，不进右端
                    }
                    break;
                }
                default: break; // 固定(1)/指定位移(2)：外力算作反力，不进右端
            }
        }
    }
}

// 装配单元刚度和右端约束修正（主从消去 / 变换法）
// 每个局部自由度展开成保留自由度的线性组合 + 常数：u_i = Σ c_{i,p} a_p + g_i
// 刚度: K(p,q) += c_{i,p} c_{j,q} ke(i,j)      —— 即 Tᵀ ke T 的散射形式
// 右端: f(p)   -= c_{i,p} ke(i,j) g_j          —— 即 -Tᵀ ke g（含指定位移与 beta）
void Assembler::AssembleStiffnessAndConstraintCorrection(Model &model) {
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            Element& element = group.GetElement(e);
            unsigned int nd = element.GetND();

            // 单元刚度
            DenseMatrix<double> ke(nd, nd);
            element.ElementStiffness(ke);

            // 预计算每个局部自由度的展开
            std::vector<DofExpansion> exp(nd);
            for (unsigned int i = 0; i < nd; i++)
                exp[i] = GetLocalDofExpansion(element, i, model);

            for (unsigned int i = 0; i < nd; i++) {
                const DofExpansion& ei = exp[i];
                for (unsigned int j = 0; j < nd; j++) {
                    const DofExpansion& ej = exp[j];
                    double kij = ke(i, j);
                    if (kij == 0.0) continue;

                    for (const auto& p : ei.terms) {
                        // 刚度：只写上三角，避免 skyline 折叠导致的重复计入
                        for (const auto& q : ej.terms) {
                            if (p.globalEqn <= q.globalEqn)
                                (*model.K)(p.globalEqn, q.globalEqn)
                                    += p.coeff * q.coeff * kij;
                        }
                        // 右端：把 j 的常数（指定位移 / beta / 指定位移master）移到右端
                        if (ej.constant != 0.0)
                            model.force[p.globalEqn - 1]
                                -= p.coeff * kij * ej.constant;
                    }
                }
            }
        }
    }
}

void Assembler::WriteDisplacementToNodes(Model &model) {
    for (auto& node : model.nodes) {
        node.UpdateNodeDisplacement(model.force);
    }
}

void Assembler::CalculateNodalBCForce(Model &model) {
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            Element& element = group.GetElement(e);
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
                const auto* c = dynamic_cast<const ContinuumElement*>(&g.GetElement(e));
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
            Element& element = group.GetElement(e);
            // 只有连续介质单元才需要外推
            auto* continuum = dynamic_cast<ContinuumElement*>(&element);
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
                Node* node = elemNodes[n];
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

// 求解从自由度位移：u_slave = Σ coeff·u_master + beta
void Assembler::RecoverSlaveDisplacement(Model &model) {
    for (auto& mpc: model.mpcs) {
        double& slaveDisp = model.nodes[mpc.slaveNode_0].displacement[mpc.slaveDof_0];
        slaveDisp += mpc.beta;
        for (auto& m: mpc.masters) {
            slaveDisp += m.coeff * model.nodes[m.node_0].displacement[m.dof_0];
        }
    }
}

// 局部自由度展开为全局自由度的线性组合+常数
// 四种 bcode 对应变换矩阵 L 的一行 + g 的一个分量
DofExpansion Assembler::GetLocalDofExpansion(const Element &element,
    unsigned int localDof, const Model &model) {
    DofExpansion result;
    // 局部自由度转化为(全局节点,全局分量)
    unsigned int ndof = element.GetNumActiveDOFsPerNode();
    const DOFIndex* dofs = element.GetActiveDOFs();
    unsigned int nodeIdx = localDof / ndof;
    unsigned int dof_0 = dofs[localDof % ndof];
    const Node* node = element.GetNodes()[nodeIdx];
    switch (node->bcode[dof_0]) {
        case 0: // 自由：展开为自身，系数为1
            result.terms.push_back({node->eqn[dof_0],1.0});
            break;
        case 1: // 固定无需展开
            break;
        case 2: // 指定位移，只有常数
            result.constant = node->displacement[dof_0];
            break;
        case 3: { // 从自由度，展开为主自由度的线性组合
            int mpcId = model.FindMPCBySlave(node->Index, dof_0);
            if (mpcId < 0) {
                throw std::runtime_error("Slave Dof has no MPC entry (node "
                    + std::to_string(node->Index+1 )+")");
            }
            const MPC& mpc = model.mpcs[mpcId];
            result.constant = mpc.beta;
            for (const auto& t: mpc.masters) {
                const Node& m = model.nodes[t.node_0];
                switch (m.bcode[t.dof_0]) {
                    case 0: // 主自由度为自由，正常进入系数中
                        result.terms.push_back({m.eqn[t.dof_0],t.coeff});
                        break;
                    case 1: // 主自由度固定，则无贡献
                        break;
                    case 2: // 主自由度指定位移，则对常数项有贡献
                        result.constant += t.coeff * m.displacement[t.dof_0];
                        break;
                    case 3: // 主自由度为为其余MPC的从自由度，则报错
                        throw std::runtime_error("Chained MPC not supported");
                    default:break;
                }
            }
            break;
        }
        default:
            throw std::runtime_error("Unknow bcode");
    }
    return result;
}

// 收集单元展开后会被写入 K 的全部全局方程号
// 无 MPC 时等价于 LocationMatrix 的非零项
std::vector<unsigned int> Assembler::GetEffectiveEquations(const Element& element,
    const Model& model) {
    std::vector<unsigned int> eqs;
    unsigned int nd = element.GetND();
    for (unsigned int i = 0; i < nd; i++) {
        DofExpansion e = GetLocalDofExpansion(element, i, model);
        for (const auto& t : e.terms)
            eqs.push_back(t.globalEqn);   // terms 里的 globalEqn 非零，否在根本不会进入循环中
    }
    return eqs;
}
//
// Created by Administrator on 2026/7/6.
//

#include "VtuExporter.h"
#include "VtuWriter.h"
#include "../Model/Model.h"
#include "../Model/Node.h"
#include "../Model/Element/Element.h"
#include "../Model/Element/ElementGroup.h"
#include "../Model/Element/CContinuumElement.h"
#include "../Core/Types.h"

namespace {
    int VtkCellType(ElementTypes t) {
        switch (t) {
            case ElementTypes::Bar3D: return 3;
            case ElementTypes::Q4_PS:
            case ElementTypes::Q4_PE: return 9;
            case ElementTypes::T3_PS:
            case ElementTypes::T3_PE:return 5;
            case ElementTypes::H8: return 12;
            case ElementTypes::Tet4: return 10;
            default: return 0;
        }
    }
    // miss应力计算函数
double VonMises(const std::vector<double>& s) {
    if (s.size() == 3) {
        double sxx = s[0], syy = s[1], sxy = s[2];
        return std::sqrt(sxx*sxx - sxx*syy + syy*syy + 3.0*sxy*sxy);
    }
    if (s.size() == 6) {
        double sxx=s[0], syy=s[1], szz=s[2];
        double sxy=s[3], syz=s[4], sxz=s[5];
        double d1 = sxx-syy, d2 = syy-szz, d3 = szz-sxx;
        return std::sqrt(0.5*(d1*d1+d2*d2+d3*d3)
                        + 3.0*(sxy*sxy+syz*syz+sxz*sxz));
    }
    return 0.0;
}

// 应力分量名称
const char* StressComponentName(unsigned int nComp, unsigned int c) {
    static const char* names2d[3] = {"Sxx", "Syy", "Sxy"};
    static const char* names3d[6] = {"Sxx", "Syy", "Szz", "Sxy", "Syz", "Sxz"};
    if (nComp == 3) return names2d[c];
    if (nComp == 6) return names3d[c];
    return "S?";
}

// 将局部坐标压平为一维度数组
std::vector<double> CollectNodePositions(const Model& model) {
    std::vector<double> xyz;
    xyz.reserve(3 * model.nodes.size());
    for (const auto& n : model.nodes) {
        xyz.push_back(n.XYZ[0]);
        xyz.push_back(n.XYZ[1]);
        xyz.push_back(n.XYZ[2]);
    }
    return xyz;
}

// 单元节点连接关系
struct Connectivity {
    std::vector<unsigned int> conn;
    std::vector<unsigned int> offs;
    std::vector<unsigned int> types;
};
Connectivity CollectConnectivity(const Model& model) {
    Connectivity out;
    int running_offset = 0;
    for (const auto& g : model.groups) {
        int vtkType = VtkCellType(g.GetElementType());
        for (unsigned int e = 0; e < g.GetNUME(); ++e) {
            const auto& elem = g.GetElement(e);
            const auto& nodes = elem.GetNodes();
            for (auto* np : nodes)
                out.conn.push_back(np->NodeNumber - 1);
            running_offset += static_cast<int>(nodes.size());
            out.offs.push_back(running_offset);
            out.types.push_back(vtkType);
        }
    }
    return out;
}

// 收集节点位移数据
std::vector<double> CollectDisplacement(const Model& model) {
    std::vector<double> u;
    u.reserve(3 * model.nodes.size());
    for (const auto& n : model.nodes) {
        u.push_back(n.Displacement[0]);
        u.push_back(n.Displacement[1]);
        u.push_back(n.Displacement[2]);
    }
    return u;
}

template <class Extractor>
std::vector<double> CollectNodeField3(const Model& model, Extractor ex) {
    std::vector<double> out;
    out.reserve(3 * model.nodes.size());
    for (const auto& n : model.nodes) {
        auto v = ex(n);   // 期望返回 std::array<double,3> 或类似
        out.push_back(v[0]);
        out.push_back(v[1]);
        out.push_back(v[2]);
    }
    return out;
}

// 节点应力
std::vector<double> CollectNodalStressComponent(const Model& model, unsigned int c) {
    std::vector<double> out;
    out.reserve(model.nodes.size());
    for (const auto& n : model.nodes) {
        out.push_back(n.stress.size() > c ? n.stress[c] : 0.0);
    }
    return out;
}
// 节点miss应力
std::vector<double> CollectNodalMises(const Model& model) {
    std::vector<double> out;
    out.reserve(model.nodes.size());
    for (const auto& n : model.nodes)
        out.push_back(VonMises(n.stress));
    return out;
}
}

bool VtuExporter::ExportMesh(const std::string& fileName, const Model& model) {
    VtuWriter w(fileName);
    if (!w.IsOpen()) return false;

    // ---- 几何 ----
    auto xyz  = CollectNodePositions(model);
    auto cell = CollectConnectivity(model);

    w.BeginPiece(model.nodes.size(),
                 cell.offs.size());
    w.WritePoints(xyz);
    w.WriteCells(cell.conn, cell.offs, cell.types);

    // ---- 节点场 ----
    w.BeginPointData();

    // Displacement (vector)
    w.WriteVectorField("Displacement", CollectDisplacement(model), 3);

    // NodeForce (vector)
    {
        std::vector<double> f;
        f.reserve(3 * model.nodes.size());
        for (const auto& n : model.nodes) {
            f.push_back(n.NodeForce[0]);
            f.push_back(n.NodeForce[1]);
            f.push_back(n.NodeForce[2]);
        }
        w.WriteVectorField("NodeForce", f, 3);
    }

    // BCForce
    {
        std::vector<double> f;
        f.reserve(3 * model.nodes.size());
        for (const auto& n : model.nodes) {
            f.push_back(n.NodeBCForce[0]);
            f.push_back(n.NodeBCForce[1]);
            f.push_back(n.NodeBCForce[2]);
        }
        w.WriteVectorField("BCForce", f, 3);
    }

    // 应力分量 + Mises（如果节点有应力）
    if (!model.nodes.empty() && !model.nodes[0].stress.empty()) {
        unsigned int nComp = model.nodes[0].stress.size();
        for (unsigned int c = 0; c < nComp; ++c) {
            w.WriteScalarField(StressComponentName(nComp, c),
                               CollectNodalStressComponent(model, c));
        }
        w.WriteScalarField("VonMises", CollectNodalMises(model));
    }

    w.EndPointData();

    // 输出单元平均miss应力，单元总应变能
    w.BeginCellData();
    // 单元编号
    {
        std::vector<unsigned int> ids;
        for (const auto& group : model.groups) {
            for (unsigned int e = 0; e < group.GetNUME(); ++e) {
                const auto& elem = group.GetElement(e);
                ids.push_back(elem.GetElementNumber());
            }
        }
        w.WriteIntField("ElementID", ids);
    }
    // 单元类型
    {
        std::vector<unsigned int> types;
        for (const auto& group : model.groups) {
            unsigned int eType = static_cast<unsigned int>(group.GetElementType());
            for (unsigned int e = 0; e < group.GetNUME(); ++e) {
                types.push_back(eType);
            }
        }
        w.WriteIntField("ElementType", types);
    }
    // 输出单元平均miss应力, Bar单元为单元内力，连续介质单元为平均miss应力
    {
        std::vector<double> s;
        for (const auto& g : model.groups)
            for (unsigned int e = 0; e < g.GetNUME(); ++e)
                s.push_back(g.GetElement(e).GetRepresentativeStress());
        w.WriteScalarField("ElementStress", s);
    }
    // 输出单元内总应变能
    {
        std::vector<double> energies;
        for (const auto& g : model.groups)
            for (unsigned int e = 0; e < g.GetNUME(); ++e)
                energies.push_back(g.GetElement(e).CalculateElementEnergy());
        w.WriteScalarField("ElementEnergy", energies);
    }
    w.EndCellData();


    w.EndPiece();

    return true;
}

bool VtuExporter::ExportGaussPoints(const std::string& fileName,
                                    const Model& model) {
    VtuWriter w(fileName);
    if (!w.IsOpen()) return false;

    // ---- 收集所有 GP ----
    std::vector<double> xyz;
    std::vector<unsigned int>    elemIDs;
    std::vector<unsigned int>    gpIDs;
    std::vector<std::vector<double>> stresses;   // 每个 GP 一条应力
    unsigned int nComp = 0;

    for (const auto& g : model.groups) {
        for (unsigned int e = 0; e < g.GetNUME(); ++e) {
            const auto& elem = g.GetElement(e);
            const auto* c    = dynamic_cast<const CContinuumElement*>(&elem);
            if (!c) continue;

            auto positions = c->GetIntegrationPointPositions();
            auto sigmas    = c->GetIntegrationPointStresses();

            const unsigned int nGp = sigmas.size();
            const unsigned int nDim = positions.GetRow();

            for (unsigned int gp = 0; gp < nGp; ++gp) {
                xyz.push_back(positions(0, gp));
                xyz.push_back(nDim >= 2 ? positions(1, gp) : 0.0);
                xyz.push_back(nDim == 3 ? positions(2, gp) : 0.0);
                elemIDs.push_back(static_cast<int>(elem.GetElementNumber()));
                gpIDs.push_back(static_cast<int>(gp + 1));
                stresses.push_back(sigmas[gp]);
                if (nComp == 0) nComp = sigmas[gp].size();
            }
        }
    }

    if (xyz.empty()) return false;

    const std::size_t nGpTotal = elemIDs.size();

    // ---- 几何 ----
    // 每个 GP 是一个 VTK_VERTEX（type=1）
    std::vector<unsigned int> conn(nGpTotal), offs(nGpTotal), types(nGpTotal, 1);
    for (unsigned int i = 0; i < nGpTotal; ++i) {
        conn[i] = i;
        offs[i] = i + 1;
    }

    w.BeginPiece(nGpTotal, nGpTotal);
    w.WritePoints(xyz);
    w.WriteCells(conn, offs, types);

    // ---- 场 ----
    w.BeginPointData();
    w.WriteIntField("ElementID",    elemIDs);
    w.WriteIntField("GaussPointID", gpIDs);

    // 应力分量
    for (unsigned int c = 0; c < nComp; ++c) {
        std::vector<double> col(nGpTotal);
        for (std::size_t i = 0; i < nGpTotal; ++i) col[i] = stresses[i][c];
        w.WriteScalarField(StressComponentName(nComp, c), col);
    }

    // Mises
    std::vector<double> mises(nGpTotal);
    for (std::size_t i = 0; i < nGpTotal; ++i) mises[i] = VonMises(stresses[i]);
    w.WriteScalarField("VonMises", mises);

    w.EndPointData();
    w.EndPiece();

    return true;
}
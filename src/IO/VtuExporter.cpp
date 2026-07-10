//
// Created by Administrator on 2026/7/6.
//

#include <fstream>
#include <iostream>
#include <cmath>
#include "VtuExporter.h"
#include "../Model/Model.h"
#include "../Model/Node.h"
#include "../Model/Element/Element.h"
#include "../Model/Element/ElementGroup.h"
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
}

bool VtuExporter::ExportVtk(const std::string &fileName, const Model &model) {
    std::ofstream out(fileName);
    if (!out.is_open()) {
        std::cerr << "Unable to open file " << fileName << std::endl;
        return false;
    }
    const unsigned int numnp = model.nodes.size();
    unsigned int totalCells = 0, totalConn = 0;
    for (auto& group : model.groups) {
        for (unsigned int e = 0; e < group.GetNUME(); e++) {
            totalCells++;
            totalConn += group.GetElement(e).GetNEN();
        }
    }
    // 写头
    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "  <UnstructuredGrid>\n"
        << "    <Piece NumberOfPoints=\"" << numnp
        << "\" NumberOfCells=\"" << totalCells << "\">\n";
    // 节点坐标
    out << "      <Points>\n"
        << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& n : model.nodes)
        out << "          " << n.XYZ[0] << " " << n.XYZ[1] << " " << n.XYZ[2] << "\n";
    out << "        </DataArray>\n      </Points>\n";

    // 单元连接
    out << "      <Cells>\n"
        << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& group : model.groups) {
        for (std::size_t e = 0; e < group.GetNUME(); ++e) {
            const auto& elem  = group.GetElement(e);
            const auto& nodes = elem.GetNodes();
            out << "         ";
            for (auto* np : nodes) out << " " << (np->NodeNumber - 1);
            out << "\n";
        }
    }

    out << "        </DataArray>\n";
    // offsets
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    std::size_t off = 0;
    for (const auto& group : model.groups) {
        for (std::size_t e = 0; e < group.GetNUME(); ++e) {
            off += group.GetElement(e).GetNEN();
            out << "          " << off << "\n";
        }
    }
    out << "\n        </DataArray>\n";
    // 单元类型
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (const auto& group : model.groups) {
        int vtkType = VtkCellType(group.GetElementType());
        for (std::size_t e = 0; e < group.GetNUME(); ++e)
            out << "          " << vtkType << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </Cells>\n";

    out << "      <PointData>\n";

    // Displacement
    out << "        <DataArray type=\"Float64\" Name=\"Displacement\" "
           "NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& n : model.nodes)
        out << "          " << n.Displacement[0] << " "
                            << n.Displacement[1] << " "
                            << n.Displacement[2] << "\n";
    out << "        </DataArray>\n";

    // NodeForce
    out << "        <DataArray type=\"Float64\" Name=\"NodeForce\" "
           "NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& n : model.nodes)
        out << "          " << n.NodeForce[0] << " "
                            << n.NodeForce[1] << " "
                            << n.NodeForce[2] << "\n";
    out << "        </DataArray>\n";

    // NodeBCForce
    out << "        <DataArray type=\"Float64\" Name=\"NodeBCForce\" "
           "NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& n : model.nodes)
        out << "          " << n.NodeBCForce[0] << " "
                            << n.NodeBCForce[1] << " "
                            << n.NodeBCForce[2] << "\n";
    out << "        </DataArray>\n";

    // NodeStress
    if (!model.nodes.empty() && !model.nodes[0].stress.empty()) {
        const unsigned int nComp = model.nodes[0].stress.size();
        const char* nodes2d[3] = {"Sxx", "Syy", "Sxy"};
        const char* nodes3d[6] = {"Sxx", "Syy" ,"Szz", "Sxy", "Syz", "Szx"};
        const char** nodesnd = (nComp == 3) ? nodes2d : nodes3d;
        for (unsigned int c = 0; c < nComp; c++) {
            out << "        <DataArray type=\"Float64\" Name=\"" << nodesnd[c]
            <<"\" format=\"ascii\">\n";
            for (auto& node : model.nodes) {
                out << "          " << (node.stress.size() > c ? node.stress[c] : 0.0) << "\n";
            }
            out << "        </DataArray>\n";
        }
    }
    // Vonmiss
    out << "        <DataArray type=\"Float64\" Name=\"VonMises\" format=\"ascii\">\n";
    for (const auto& n : model.nodes) {
        double vm = 0.0;
        if (n.stress.size() == 3) {
            // 2D 平面应力，后续添加平面应变状态
            double sxx = n.stress[0], syy = n.stress[1], sxy = n.stress[2];
            vm = std::sqrt(sxx*sxx - sxx*syy + syy*syy + 3.0*sxy*sxy);
        } else if (n.stress.size() == 6) {
            // 3D
            double sxx = n.stress[0], syy = n.stress[1], szz = n.stress[2];
            double sxy = n.stress[3], syz = n.stress[4], sxz = n.stress[5];
            double s1 = sxx - syy, s2 = syy - szz, s3 = szz - sxx;
            vm = std::sqrt(0.5*(s1*s1 + s2*s2 + s3*s3) + 3.0*(sxy*sxy + syz*syz + sxz*sxz));
        }
        out << "          " << vm << "\n";
    }
    out << "        </DataArray>\n";

    out << "      </PointData>\n";
    out << "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
    return true;
}

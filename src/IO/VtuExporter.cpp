//
// Created by Administrator on 2026/7/6.
//

#include "VtuExporter.h"
#include "../Model/Model.h"
#include "../Model/Node.h"
#include "../Model/Element/Element.h"
#include "../Model/Element/ElementGroup.h"
#include "../Core/Types.h"
#include <fstream>
#include <iostream>

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
            out << " " << off;
        }
    }
    out << "\n        </DataArray>\n";
    // 单元类型
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (const auto& group : model.groups) {
        int vtkType = VtkCellType(group.GetElementType());
        for (std::size_t e = 0; e < group.GetNUME(); ++e)
            out << " " << vtkType;
    }
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

    out << "      </PointData>\n";
    out << "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
    return true;
}

//
// Created by Administrator on 2026/7/10.
//

#include "VtuWriter.h"

#include <iostream>

VtuWriter::VtuWriter(const std::string &fileName): out_(fileName){
    if (!out_.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return;
    }
    // 写头
    out_ << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "  <UnstructuredGrid>\n";
}

VtuWriter::~VtuWriter() {
    if (out_.is_open()) {
        out_ << "  </UnstructuredGrid>\n</VTKFile>\n";
    }
}

bool VtuWriter::IsOpen() {
    if (out_.is_open()) {
        return true;
    }
    return false;
}

void VtuWriter::BeginPiece(unsigned int nP, unsigned int nC) {
    out_ << "    <Piece NumberOfPoints=\"" << nP
         << "\" NumberOfCells=\"" << nC << "\">\n";
}
void VtuWriter::EndPiece() {
    out_ << "    </Piece>\n";
}

void VtuWriter::WritePoints(const std::vector<double>& xyz) {
    out_ << "      <Points>\n"
         << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < xyz.size(); i += 3) {
        out_ << "          "
             << xyz[i]   << " "
             << xyz[i+1] << " "
             << xyz[i+2] << "\n";
    }
    out_ << "        </DataArray>\n"
         << "      </Points>\n";
}

void VtuWriter::WriteCells(const std::vector<unsigned int>& conn,
                           const std::vector<unsigned int>& offs,
                           const std::vector<unsigned int>& types) {
    out_ << "      <Cells>\n";

    out_ << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int v : conn) out_ << "          " << v << "\n";
    out_ << "        </DataArray>\n";

    out_ << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int v : offs) out_ << "          " << v << "\n";
    out_ << "        </DataArray>\n";

    out_ << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int v : types) out_ << "          " << v << "\n";
    out_ << "        </DataArray>\n";

    out_ << "      </Cells>\n";
}

void VtuWriter::BeginPointData() { out_ << "      <PointData>\n"; }
void VtuWriter::EndPointData()   { out_ << "      </PointData>\n"; }
void VtuWriter::BeginCellData()  { out_ << "      <CellData>\n"; }
void VtuWriter::EndCellData()    { out_ << "      </CellData>\n"; }

void VtuWriter::WriteScalarField(const std::string& name,
                                 const std::vector<double>& data) {
    out_ << "        <DataArray type=\"Float64\" Name=\"" << name
         << "\" format=\"ascii\">\n";
    for (double v : data) out_ << "          " << v << "\n";
    out_ << "        </DataArray>\n";
}

void VtuWriter::WriteVectorField(const std::string& name,
                                 const std::vector<double>& data,
                                 unsigned int nComp) {
    out_ << "        <DataArray type=\"Float64\" Name=\"" << name
         << "\" NumberOfComponents=\"" << nComp
         << "\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < data.size(); i += nComp) {
        out_ << "          ";
        for (unsigned int c = 0; c < nComp; ++c)
            out_ << data[i + c] << " ";
        out_ << "\n";
    }
    out_ << "        </DataArray>\n";
}

void VtuWriter::WriteIntField(const std::string& name,
                              const std::vector<unsigned int>& data) {
    out_ << "        <DataArray type=\"Int32\" Name=\"" << name
         << "\" format=\"ascii\">\n";
    for (int v : data) out_ << "          " << v << "\n";
    out_ << "        </DataArray>\n";
}
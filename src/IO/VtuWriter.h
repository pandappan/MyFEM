#include <fstream>
#include <string>
#include <vector>

#pragma once
class VtuWriter {
private:
    std::ofstream out_;
    int indent_;

public:
    explicit VtuWriter(const std::string& fileName);
    ~VtuWriter();
    bool IsOpen();
    void BeginPiece(unsigned int numPoints, unsigned int numCells);
    void EndPiece();

    void WritePoints(const std::vector<double>& xyz);
    void WriteCells(
        const std::vector<unsigned int> &conn,
        const std::vector<unsigned int> &offs,
        const std::vector<unsigned int> &types);

    void BeginPointData();
    void EndPointData();

    void BeginCellData();
    void EndCellData();

    void WriteScalarField(const std::string& name, const std::vector<double>& data);
    void WriteVectorField(const std::string& name,
        const std::vector<double>& data,
        unsigned int numComponents);
    void WriteIntField(const std::string &name, const std::vector<unsigned int> &data);
};

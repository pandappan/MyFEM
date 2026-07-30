#include <iomanip>
#include "AxisymMaterial.h"
#include "Types.h"
#include "DenseMatrix.h"

AxisymMaterial::AxisymMaterial() {
    matType = MaterialCategory::MechanicalAxisym;
}

void AxisymMaterial::Write(std::ostream &output) const {
    output << std::setw(6) << rho
    << std::setw(6) << E
    << std::setw(6) << nu << std::endl;
}

// 应变 {εrr, εzz, εθθ, γrz}, 应力 {σrr, σzz, σθθ, τrz}
void AxisymMaterial::ComputeElasticMatrix(DenseMatrix<double>& D) const {
    D.SetZero();
    double a = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double b = E * nu         / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double G = E / (2.0 * (1.0 + nu));

    D(0,0) = a; D(0,1) = b; D(0,2) = b;
    D(1,0) = b; D(1,1) = a; D(1,2) = b;
    D(2,0) = b; D(2,1) = b; D(2,2) = a;
    D(3,3) = G;
}


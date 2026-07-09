//
// Created by Administrator on 2026/7/6.
//

#include "Writer.h"
#include "../Model/Model.h"
#include "../Model/Node.h"
#include "../Model/Element/Element.h"
#include "../Model/Element/ElementGroup.h"
#include "../Model/Material/Material.h"
#include "../Core/SkylineMatrix.h"
#include <ctime>
#include <cstdlib>

Writer::Writer(const std::string& fileName) : outputFile_(fileName) {
    if (!outputFile_.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        std::exit(3);
    }
}

void Writer::PrintTime() {
    time_t rawtime;
    time(&rawtime);
    tm* t = localtime(&rawtime);
    static const char* wd[] = {"Sunday","Monday","Tuesday","Wednesday",
                               "Thursday","Friday","Saturday"};
    static const char* mo[] = {"January","February","March","April","May","June",
                               "July","August","September","October","November","December"};
    *this << "        ("
          << t->tm_hour << ":" << t->tm_min << ":" << t->tm_sec
          << " on " << mo[t->tm_mon] << " " << t->tm_mday << ", "
          << t->tm_year + 1900 << ", " << wd[t->tm_wday] << ")\n\n";
}

void Writer::OutputHeading(const Model &model) {
    *this << "Title: " << model.title << std::endl;
    PrintTime();
}

void Writer::OutputNodeInfo(const Model& model) {
    *this << "C O N T R O L   I N F O R M A T I O N\n\n";
    *this << "  DIMENSION = " << model.dimension << "\n";
    *this << "  NUMNP     = " << model.nodes.size() << "\n";
    *this << "  NUMEG     = " << model.groups.size() << "\n\n";
    Tee([&](std::ostream& os) {
        for (auto& node : model.nodes) node.Write(os, model.dimension);
    });
    *this << "\n";
}

void Writer::OutputEquationNumber(const Model& model) {
    *this << " EQUATION NUMBERS\n\n"
          << "   NODE   Degrees of freedom\n";
    if (model.dimension == 2)
        *this << "     N     X    Y\n";
    else
        *this << "     N     X    Y    Z\n";
    Tee([&](std::ostream& os) {
        for (auto& n : model.nodes) n.WriteEquationNo(os, model.dimension);
    });
    *this << "\n";
}
void Writer::OutputElementInfo(const Model& model) {
    *this << " E L E M E N T   G R O U P   D A T A\n\n\n";
    for (const auto & group : model.groups) {
        *this << " E L E M E N T   D E F I N I T I O N\n\n";
        *this << "   ELEMENT TYPE  = " << static_cast<int>(group.GetElementType()) << "\n";
        *this << "   NUME          = " << group.GetNUME()   << "\n";
        *this << "   NUMMAT        = " << group.GetNUMAT()  << "\n\n";
        // 材料列表
        *this << " M A T E R I A L   D E F I N I T I O N\n";
        Tee([&](std::ostream& os) {
            for (std::size_t m = 0; m < group.GetNUMAT(); ++m) {
                os << std::setw(5) << m + 1;
                group.GetMaterial(m).Write(os);
            }
        });
        // 单元列表
        *this << "\n E L E M E N T   I N F O R M A T I O N\n";
        Tee([&](std::ostream& os) {
            for (std::size_t e = 0; e < group.GetNUME(); ++e) {
                os << std::setw(5) << e + 1;
                group.GetElement(e).Write(os);
            }
        });
        *this << "\n";
    }
}

void Writer::OutputNodeForce(const Model& model) {
    *this << " N O D A L   F O R C E S\n\n"
          << "  NODE      ";
    if (model.dimension == 2)
        *this << "        FX                FY\n";
    else
        *this << "        FX                FY                FZ\n";
    Tee([&](std::ostream& os) {
        for (auto& n : model.nodes) n.WriteNodeForces(os, model.dimension);
    });
    *this << "\n";
}

void Writer::OutputNodalDisplacement(const Model& model) {
    *this << " D I S P L A C E M E N T S\n\n";
    if (model.dimension == 2)
        *this << "  NODE      X-DISPLACEMENT    Y-DISPLACEMENT\n";
    else
        *this << "  NODE      X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT\n";
    Tee([&](std::ostream& os) {
        for (auto& n : model.nodes) n.WriteNodalDisplacement(os, model.dimension);
    });
    *this << "\n";
}

void Writer::OutputNodalBCForce(const Model &model) {
    *this << "BCForce\n\n";
    if (model.dimension == 2) {
        *this << " NODE    X-BCForce    Y-BCForce\n";
    } else {
        *this << " NODE    X-BCForce    Y-BCForce    Z-BCForce\n";
    }
        Tee([&](std::ostream& os) {
        for (auto& n : model.nodes) n.WriteNodeBCForces(os, model.dimension);
    });
    *this << "\n";
}

void Writer::OutputElementStress(const Model& model) {
    unsigned int i = 0;
    for (const auto& group : model.groups) {
        *this << " STRESSES FOR ELEMENT GROUP " << i++ << "\n\n";
        Tee([&](std::ostream& os) {
            for (std::size_t e = 0; e < group.GetNUME(); ++e) {
                group.GetElement(e).WriteElementStress(os);
            }
        });
        *this << "\n";
    }
}

void Writer::OutputTotalSystemData(const Model& model) {
    if (!model.K) return;
    *this << "\nTOTAL SYSTEM DATA\n"
          << "  NEQ = " << model.neq << "\n"
          << "  NWK = " << model.K->size() << "\n"
          << "  MK  = " << model.K->GetMaximumHalfBandwidth() << "\n";
}

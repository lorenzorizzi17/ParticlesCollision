#include "ParticleType.hpp"



ParticleType::ParticleType(std::string const& Name, double Mass, int Charge) : fName{Name}, fMass{Mass}, fCharge{Charge} {}
std::string ParticleType::GetName() const {
    return fName;
};
double ParticleType::GetMass() const {
    return fMass;
};
int ParticleType::GetCharge() const {
    return fCharge;
};

void ParticleType::Print() const {
    std::cout << "Particle " << fName << " of mass " << fMass << ", charge " << fCharge;
};

double ParticleType::GetWidth() const {
    return 0;
};

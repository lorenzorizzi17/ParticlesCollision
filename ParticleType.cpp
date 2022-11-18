#include "ParticleType.hpp"


//parametric constructor
ParticleType::ParticleType(std::string const& Name, double Mass, int Charge) : fName{Name}, fMass{Mass}, fCharge{Charge} {}

//getters
std::string ParticleType::GetName() const {
    return fName;
};
double ParticleType::GetMass() const {
    return fMass;
};
int ParticleType::GetCharge() const {
    return fCharge;
};

//print method 
void ParticleType::Print() const {
    std::cout << "Particle " << fName << " of mass " << fMass << ", charge " << fCharge;
};

//width getter (if not overridden, it returns 0)
double ParticleType::GetWidth() const {
    return 0;
};

#include"ResonanceType.hpp"

//parametric constructor
ResonanceType::ResonanceType(std::string const& Name, double Mass, int Charge, double Width): ParticleType(Name, Mass, Charge), fWidth{Width} {};

//width getter
double ResonanceType::GetWidth() const {return fWidth;}

//print method overriden
void ResonanceType::Print() const {
    ParticleType::Print();
    std::cout << " and width " << fWidth;
}

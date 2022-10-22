#include"ResonanceType.hpp"

ResonanceType::ResonanceType(std::string const& Name, double Mass, int Charge, double Width): ParticleType(Name, Mass, Charge), fWidth{Width} {};
double ResonanceType::GetWidth() const {return fWidth;}
void ResonanceType::Print() const {
    ParticleType::Print();
    std::cout << " and width " << fWidth;
}

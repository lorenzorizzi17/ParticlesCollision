#ifndef RESONANCE_TYPE
#define RESONANCE_TYPE
#include "ParticleType.hpp"

class ResonanceType : public ParticleType {
    public:
    ResonanceType(std::string const&, double, int, double);
    double GetWidth() const override;
    void Print() const override;
    private:
    double const fWidth;

};

#endif
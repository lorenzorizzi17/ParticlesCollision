#ifndef PARTICLE_TYPE
#define PARTICLE_TYPE

#include <string>
#include <iostream>

class ParticleType {
    public:
    ParticleType(std::string const&, double, int);
    std::string GetName() const;
    double GetMass() const;
    int GetCharge() const;
    virtual void Print() const;
    virtual double GetWidth() const;


    private:
    const std::string fName;
    const double fMass;
    const int fCharge;
};

#endif
#ifndef PARTICLE_H
#define PARTICLE_H

#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include <vector>
#include <cmath>

class Particle {
    public:
    Particle(std::string const&, double d = 0, double e = 0, double f = 0);
    Particle();
    static void AddParticleType(std::string const&, double, int, double);
    int GetIndex() const;
    int GetCharge() const;
    std::string GetName() const;
    void PrintParticleInfo() const;
    void SetIndex(int);
    void SetIndex(std::string const&);
    static void PrintParticleVector();
    double GetPx() const;
    double GetPy() const;
    double GetPz() const;
    double GetMass() const;
    double GetEnergy() const;
    double InvMass(Particle const& p) const;
    void SetP(double, double, double);
    int Decay2body(Particle &dau1,Particle &dau2) const;
  


    private:
    static std::vector<ParticleType*> fParticleType;
    static int fMaxNumParticleType;
    int fIndex;
    double fPx, fPy, fPz;
    static int FindParticle(std::string const&);
    void Boost(double bx, double by, double bz);
};


#endif
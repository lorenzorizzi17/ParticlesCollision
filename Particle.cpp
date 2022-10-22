#include "Particle.hpp"
#include <iterator>
#include <algorithm>
#include <math.h>  

std::vector<ParticleType *> Particle::fParticleType{};

Particle::Particle(std::string const &Name, double Px, double Py, double Pz) : fPx{Px}, fPy{Py}, fPz{Px}
{
    fIndex = FindParticle(Name);
}

int Particle::FindParticle(std::string const &s)
{
    // algorithm
    int counter = 0;
    for (auto it = fParticleType.begin(); it != fParticleType.end(); ++it)
    {
        if (s == (*it)->GetName())
        {
            return counter;
        }
        counter++;
    }
    throw std::runtime_error{"The desired particle does not exist. Try again"};
};

void Particle::AddParticleType(std::string const &Name, double Mass, int Charge, double Width)
{
    auto it = std::find_if(fParticleType.begin(), fParticleType.end(), [&](ParticleType *p)
                           { return p->GetName() == Name; });
    if (it == fParticleType.end())
    {
        if (Width > 0)
        {
            ResonanceType *resonanceType = new ResonanceType(Name, Mass, Charge, Width);
            fParticleType.push_back(resonanceType);
        }
        else
        {
            ParticleType *particleType = new ParticleType(Name, Mass, Charge);
            fParticleType.push_back(particleType);
        };
    } else {
        std::cout << "\nParticle Type already exists\n";
    }
}

int Particle::GetIndex() const
{
    return fIndex;
}

void Particle::SetIndex(int Index)
{
    
    if (Index <= fParticleType.size())
    {
        fIndex = Index;
    }
    else
    {
        std::cout << "\nParticles does not exist\n";
    }
}

void Particle::SetIndex(std::string const &s)
{
    auto it = std::find_if(fParticleType.begin(), fParticleType.end(), [&](ParticleType *p)
                           { return p->GetName() == s; });
    if (it == fParticleType.end()) {
        std::cout << "\nThere's no existing " << s << " particle\n";
    } else {
        int n = FindParticle(s);
        fIndex = n;
    }
}

void Particle::PrintParticleVector()
{
    std::cout << "\nYou requested the particle type list\n";
    for (auto it = fParticleType.begin(); it != fParticleType.end(); ++it)
    {
        (*it)->Print();
        std::cout << '\n';
    }
}

void Particle::PrintParticleInfo() const
{
    std::cout << "*******************************************************" << '\n';
    std::cout << "Particle index: " << fIndex << '\n'
              << "Particle's name: " << fParticleType[fIndex]->GetName();
    std::cout << '\n'
              << "Particle's momentum: (" << fPx << ", " << fPy << ", " << fPz << ")" << '\n';
    std::cout << "*******************************************************" << '\n';
}

double Particle::GetPx() const {
    return fPx;
}
double Particle::GetPy() const {
    return fPy;
}
double Particle::GetPz() const {
    return fPz;
}
std::string Particle::GetName() const
{
    return fParticleType[fIndex]->GetName();
}

double Particle::GetMass() const {
    return fParticleType[fIndex]->GetMass();
}

double Particle::GetEnergy() const {
    double modp2 = fPx*fPx + fPy*fPy + fPz*fPz;
    return sqrt(Particle::GetMass()*Particle::GetMass()+modp2);
}

double Particle::InvMass(Particle const& p) const {
    double E1 = this->GetEnergy();
    double E2 = p.GetEnergy();
    double E = (E1+E2)*(E1+E2);
    double modp2 = (fPx+p.GetPx())*(fPx+p.GetPx()) + (fPy+p.GetPy())*(fPy+p.GetPy()) + (fPz+p.GetPz())*(fPz+p.GetPz());
    return sqrt(E-modp2);
}

void Particle::SetP(double px, double py, double pz) {
        fPx = px;
        fPy = py;
        fPz = pz;
}
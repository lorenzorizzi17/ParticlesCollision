#include "Particle.hpp"
#include <iterator>
#include <algorithm>
#include <math.h>  
#include <cmath>  // for M_PI
#include <cstdlib> //for RAND_MAX

std::vector<ParticleType *> Particle::fParticleType{};

Particle::Particle() {
    fIndex = -1;
}

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

int Particle::GetCharge() const {
    return fParticleType[fIndex]->GetCharge();
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
int Particle::Decay2body(Particle &dau1,Particle &dau2) const {
  if(GetMass() == 0.0){
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }
  
  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if(fIndex > -1){ // add width effect

    // gaussian random numbers

    float x1, x2, w, y1, y2;
    
    double invnum = 1./RAND_MAX;
    do {
      x1 = 2.0 * rand()*invnum - 1.0;
      x2 = 2.0 * rand()*invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;

    massMot += fParticleType[fIndex]->GetWidth() * y1;

  }

  if(massMot < massDau1 + massDau2){
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }
  
  double pout = sqrt((massMot*massMot - (massDau1+massDau2)*(massDau1+massDau2))*(massMot*massMot - (massDau1-massDau2)*(massDau1-massDau2)))/massMot*0.5;

  double norm = 2*M_PI/RAND_MAX;

  double phi = rand()*norm;
  double theta = rand()*norm*0.5 - M_PI/2.;
  dau1.SetP(pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta));
  dau2.SetP(-pout*sin(theta)*cos(phi),-pout*sin(theta)*sin(phi),-pout*cos(theta));

  double energy = sqrt(fPx*fPx + fPy*fPy + fPz*fPz + massMot*massMot);

  double bx = fPx/energy;
  double by = fPy/energy;
  double bz = fPz/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);

  return 0;
}
void Particle::Boost(double bx, double by, double bz)
{

  double energy = GetEnergy();

  //Boost this Lorentz vector
  double b2 = bx*bx + by*by + bz*bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx*fPx + by*fPy + bz*fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

  fPx += gamma2*bp*bx + gamma*bx*energy;
  fPy += gamma2*bp*by + gamma*by*energy;
  fPz += gamma2*bp*bz + gamma*bz*energy;
}
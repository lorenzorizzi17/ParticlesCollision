#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "Particle.hpp"


TEST_CASE("Checking for ParticleType and ResonanceType") {
    ResonanceType resonanceType = ResonanceType("P+",2.5,2,1.33);
    ParticleType particleType = ParticleType("P+",3.,2);
    //creating the array
    ParticleType* arr[2]{&particleType, &resonanceType};
    // p points to particletype 
    ParticleType* p = arr[0];
    CHECK(p->GetWidth() == doctest::Approx(0));
    CHECK(p->GetCharge() == 2);
    CHECK(p->GetMass() == doctest::Approx(3));
    CHECK(p->GetName() == "P+");
    //p->Print(); std::cout << '\n';
    //p now points to resonancetype
    p = arr[1];
    CHECK(p->GetWidth() == doctest::Approx(1.33));
    CHECK(p->GetCharge() == 2);
    CHECK(p->GetMass() == doctest::Approx(2.5));
    CHECK(p->GetName() == "P+");
    //p->Print();std::cout<<'\n';
}


TEST_CASE("Checking for Particle"){
    //filling the vector of particles
    Particle::AddParticleType("P+",20.4, 1, 0);
    Particle::AddParticleType("e-", 2, -1,0);
    Particle::AddParticleType("K*", 20, 0, 12);
    //instantiating some stable particles
    Particle proton = Particle("P+");
    Particle electron1 = Particle("e-");
    //instantiating an unstable particle
    Particle kaon = Particle("K*");
    //Checking GetIndex()
    CHECK(proton.GetIndex() == 0);
    CHECK(electron1.GetIndex() == 1);
    CHECK(kaon.GetIndex() == 2);
    //Printing info about particles
    proton.PrintParticleInfo();
    electron1.PrintParticleInfo();
    kaon.PrintParticleInfo();

    //shifting particles
    proton.SetIndex("K*");
    //this shoul print a K* particle:
    proton.PrintParticleInfo();
};

TEST_CASE("Paolino"){
    //filling the vector of particles
    Particle::AddParticleType("P+",20.4, 1, 0);
    Particle::AddParticleType("e-", 2, -1,0);
    Particle::AddParticleType("K*", 20, 0, 12);
    //instantiating some stable particles
    Particle proton = Particle("P+");
    Particle electron1 = Particle("e-");
    //instantiating an unstable particle
    Particle kaon = Particle("K*");

    Particle::PrintParticleVector();
    CHECK(proton.GetPx() == 0);
    CHECK(kaon.GetPx() == 0);
    CHECK(proton.GetMass() == doctest::Approx(20.4));
    CHECK(electron1.GetMass() == doctest::Approx(2));
    CHECK(proton.GetEnergy() == doctest::Approx(proton.GetMass()));
    CHECK(kaon.GetEnergy() == doctest::Approx(kaon.GetMass()));
    CHECK(kaon.InvMass(proton) == doctest::Approx(40.4));

    proton.SetP(10,10,10);
    electron1.SetP(1,2,1);
    CHECK(proton.GetPx() == 10);
    CHECK(electron1.GetPx() == 1);
    CHECK(proton.GetPy() == 10);
    CHECK(electron1.GetPy() == 2);
    CHECK(proton.GetPz() == 10);
    CHECK(electron1.GetPz() == 1);
    CHECK(proton.GetEnergy() == doctest::Approx(26.761));
    CHECK(electron1.GetEnergy() == doctest::Approx(3.162).epsilon(0.001));
    CHECK(proton.InvMass(electron1) == doctest::Approx(22.5696).epsilon(0.0001));
} 
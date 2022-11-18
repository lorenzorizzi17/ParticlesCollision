#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "Particle.hpp"


TEST_CASE("Testing ParticleType") {
    std::cout << "\nTESTCASE 1: Testing ParticleType\n";
    //creating a ParticleType object
    ParticleType particleType = ParticleType("P+",3.,2);
    //testing ParticleType's getters
    CHECK(particleType.GetName() == "P+");
    CHECK(particleType.GetMass() == 3.);
    CHECK(particleType.GetCharge()==2);
    CHECK(particleType.GetWidth() == 0);
    //printing info about the particletype
    std::cout << "You invoked .Print() on ParticleType: \n";
    particleType.Print(); std::cout << '\n';
}


TEST_CASE("Testing ResonaceType"){
    std::cout << "\nTESTCASE 2: Testing ResonanceType\n";
    //creating a ResonaceType object
    ResonanceType resonanceType = ResonanceType("K*",4.4,3,1.25);
    //testing ResonanceType's getters
    CHECK(resonanceType.GetName() == "K*");
    CHECK(resonanceType.GetMass() == 4.4);
    CHECK(resonanceType.GetCharge()==3);
    CHECK(resonanceType.GetWidth() == 1.25);
    //printing info about the resonancetype
    std::cout << "You invoked .Print() on ResonanceType: \n";
    resonanceType.Print();std::cout << '\n';
}

TEST_CASE("Testing ParticleType and ResonanceType together for inheritance"){
    std::cout << "\nTESTCASE 3: Testing ParticleType and ResonanceType together for inheritance\n";
    //creating a ParticleType and a ResonanceType object
    ParticleType particleType = ParticleType("Pi-",4.2,1);
    ResonanceType resonanceType = ResonanceType("K*",2.5,2,1.33);
    //creating an array of pointers
    ParticleType* arr[2]{&particleType, &resonanceType};
    //creating a pointer p pointing to particletype 
    ParticleType* p = arr[0];
    CHECK(p->GetWidth() == doctest::Approx(0));
    CHECK(p->GetCharge() == 1);
    CHECK(p->GetMass() == doctest::Approx(4.2));
    CHECK(p->GetName() == "Pi-");
    //calling Print 
    std::cout << "You invoked ->Print() on a pointer to ParticleType: \n";
    p->Print(); std::cout << '\n';
    //p now points to resonancetype
    p = arr[1];
    CHECK(p->GetWidth() == doctest::Approx(1.33));
    CHECK(p->GetCharge() == 2);
    CHECK(p->GetMass() == doctest::Approx(2.5));
    CHECK(p->GetName() == "K*");
    std::cout << "You invoked ->Print() on a pointer to ResonanceType: \n";
    p->Print(); std::cout << '\n';
}


TEST_CASE("Testing Particle"){
    std::cout << "\nTESTCASE 4: Testing Particle\n";
    //creating a Particle instance with default constructor
    Particle particle1 = Particle();
    CHECK(particle1.GetIndex()==-1);
    //adding some particles to the ParticlesVector
    Particle::AddParticleType("P+",20.4, 1, 0);
    Particle::AddParticleType("e-", 2, -1,0);
    Particle::AddParticleType("K*", 20, 0, 12);
    //instantiating some stables particles
    Particle proton = Particle("P+");
    Particle electron1 = Particle("e-");
    //instantiating an unstable particle
    Particle kaon = Particle("K*");
    //testing getters
    CHECK(proton.GetIndex() == 0);
    CHECK(electron1.GetIndex() == 1);
    CHECK(kaon.GetIndex() == 2);
    CHECK(proton.GetMass() == doctest::Approx(20.4));
    CHECK(proton.GetCharge() == 1);
    CHECK(proton.GetName() == "P+");
    CHECK(proton.GetPx() == 0);
    //testing px setter
    proton.SetP(1,2,3);
    CHECK(proton.GetPx() == 1);
    CHECK(proton.GetPy() == 2);
    CHECK(proton.GetPz() == 3);
    //Printing info about particles
    std::cout << "You called PrintParticleInfo() on proton: \n";
    proton.PrintParticleInfo();
    std::cout << "You called PrintParticleInfo() on electron: \n";
    electron1.PrintParticleInfo();
    std::cout << "You called PrintParticleInfo() on kaon: \n";
    kaon.PrintParticleInfo();
    //shifting particles
    proton.SetIndex("K*");
    //this should now print a K* particle:
    std::cout << "You called PrintParticleInfo() on kaon: \n";
    proton.PrintParticleInfo();
    //testing energy getters
    proton.SetIndex("P+");
    CHECK(proton.GetEnergy() == doctest::Approx(20.74).epsilon(0.01));
    //testing invmass
    CHECK(kaon.InvMass(proton) == doctest::Approx(40.57).epsilon(0.01));
    //printing particles vector
    std::cout << "You requested the particles vector: \n";
    Particle::PrintParticleVector();
};

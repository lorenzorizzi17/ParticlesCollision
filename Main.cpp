#include"Particle.hpp"
#include"TRandom.h"
#include"TH1F.h"
#include<vector>
#include"TCanvas.h"
#include"TFile.h"
#include"TBenchmark.h"


R__LOAD_LIBRARY(ParticleType_cpp.so)
R__LOAD_LIBRARY(ResonanceType_cpp.so)
R__LOAD_LIBRARY(Particle_cpp.so)


void generate() {
    //Adding particle types
    Particle::AddParticleType("Pi+",0.13957,1,0);
    Particle::AddParticleType("Pi-",0.13957,-1,0);
    Particle::AddParticleType("K+", 0.49367,1,0);
    Particle::AddParticleType("K-", 0.49367,-1,0);
    Particle::AddParticleType("P+",0.93827,1,0);
    Particle::AddParticleType("P-",0.93827,-1,0);
    Particle::AddParticleType("K*",0.89166,0,0.050);
    //Setting the seed
    gRandom->SetSeed();
    //instantiating root objects such as file, graphic canvas and histograms
    TFile* File = new TFile("Particles.root","RECREATE");

    TCanvas* Canvas = new TCanvas("Particles.root","Particles.root", 1400,800);
    Canvas->Divide(4,2);

    TH1F* HistoPhi = new TH1F("Azimuthal angle", "Azimuthal angle", 1000,0,2*M_PI);
    TH1F* HistoTheta = new TH1F("Polar angle", "Polar angle", 1000,0,M_PI);
    TH1F* HistoParticleType = new TH1F("Particle type", "Particle type",7,0,7);
    TH1F* HistoMomentum = new TH1F("Momentum distribution", "Momentum distribution", 1000,0,6);
    TH1F* HistoTransverseMomentum = new TH1F("Transverse momentum distribution", "Transverse momentum distribution",1000,0,4);
    TH1F* HistoEnergy = new TH1F("Energy distribution", "Energy distribution", 1000, 0, 3);
    TH1F* HistoInvMass = new TH1F("Invariant mass distribution", "Invariant mass distribution", 1000,0,6);
    TH1F* HistoInvMassConcordi = new TH1F("Inv. mass distr. concordant particles", "Inv. mass distr. concordant particles",150,0,3);
    TH1F* HistoInvMassDiscordi = new TH1F("Inv. mass distr. discordant particles", "Inv. mass distr. discordant particles",150,0,3);
    TH1F* HistoInvMassPionKaonConcordant = new TH1F("Inv. mass distr. concordant Pi and K particles", "Inv. mass distr. concordant Pi and K particles", 150,0,3);
    TH1F* HistoInvMassPionKaonDiscordant = new TH1F("Inv. mass distr. discordant Pi and K particles", "Inv. mass distr. discordant Pi and K particles", 150,0,3);
    TH1F* HistoInvMassDecayed = new TH1F("Inv. mass distr. daughter particles", "Inv. mass distr. daughter particles",150,0.5,3);
    HistoInvMassConcordi->Sumw2();
    HistoInvMassDiscordi->Sumw2();
    HistoInvMassPionKaonConcordant->Sumw2();
    HistoInvMassPionKaonDiscordant->Sumw2();
    HistoInvMassDecayed->Sumw2();
    
    //particle generation beginning
    gBenchmark->Start("");
    std::array<Particle,120> Particles;
    for (int i = 0; i < 1E5; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            //generating a generic particle
            Particle particle = Particle();
            //generating and setting the particle momentum
            double phi = gRandom->Uniform(0,2*M_PI);
            double theta = gRandom->Uniform(0,M_PI);
            double pmod = gRandom->Exp(1);
            double px = pmod * std::sin(theta)*std::cos(phi);
            double py = pmod * std::sin(theta)*std::sin(phi);
            double pz = pmod * std::cos(theta); 
            particle.SetP(px,py,pz);
            //setting the particle type
            double test = gRandom->Rndm();
            if (test < 0.4) {
                particle.SetIndex("Pi+");
                HistoParticleType->Fill(0);
            } else if (test<0.8) {
                particle.SetIndex("Pi-");
                HistoParticleType->Fill(1);
            } else if (test < 0.85) {
                particle.SetIndex("K+");
                HistoParticleType->Fill(2);
            } else if (test < 0.9) {
                particle.SetIndex("K-");
                HistoParticleType->Fill(3);
            } else if (test <0.945) {
                particle.SetIndex("P+");
                HistoParticleType->Fill(4);
            } else if (test < 0.99) {
                particle.SetIndex("P-");
                HistoParticleType->Fill(5);
            } else {
                particle.SetIndex("K*");
                HistoParticleType->Fill(6);
            }
            HistoPhi->Fill(phi);
            HistoTheta->Fill(theta);
            HistoMomentum->Fill(pmod);
            HistoTransverseMomentum->Fill(std::sqrt(px*px+py*py));
            HistoEnergy->Fill(particle.GetEnergy());
            Particles[j] = particle;
        }
        
        int k = 100;           //k represents the current number of particles generated in an event
        //checking for K* particles to make them decay 
        for (int i{0}; i < 100; i++) {
            if (Particles[i].GetName() == "K*") {
                //making the K* particle decay
                double test2 = gRandom->Rndm();
                if (test2 < 0.5) {
                    Particle PioneP = Particle("Pi+");
                    Particle KaoneN = Particle("K-");
                    Particles[i].Decay2body(PioneP, KaoneN);
                    Particles[k] = PioneP;
                    Particles[k+1] = KaoneN; 
                    k = k+2;
                    HistoInvMassDecayed->Fill(PioneP.InvMass(KaoneN));
                } else {
                    Particle PioneN = Particle("Pi-");
                    Particle KaoneP = Particle("K+");
                    Particles[i].Decay2body(PioneN, KaoneP);
                    Particles[k] = PioneN;
                    Particles[k+1] = KaoneP; 
                    k = k+2;
                    HistoInvMassDecayed->Fill(PioneN.InvMass(KaoneP));
                }
            }
        }

        //computing the invariant mass ignoring K* particles that decayed
        for (int i = 0; i < k; i++)
        {
            for (int j = i; j < k; j++)
            {
                if (j!=i) {
                    double inv_mass = Particles[i].InvMass(Particles[j]);
                    if ((Particles[j].GetName() != "K*")&&(Particles[i].GetName() != "K*")) {
                        HistoInvMass->Fill(inv_mass);
                    }
                    //Here goes concordant/discordant particles and pi/k combination
                    if ((Particles[i].GetCharge()*Particles[j].GetCharge() > 0)&&(((Particles[j].GetName() != "K*")&&(Particles[i].GetName() != "K*")))){
                        HistoInvMassConcordi->Fill(inv_mass);
                        if(((Particles[i].GetName() == "Pi+")&&(Particles[j].GetName()=="K+"))||((Particles[i].GetName() == "K+")&&(Particles[j].GetName()=="Pi+"))){
                        HistoInvMassPionKaonConcordant->Fill(inv_mass);
                        }
                        if(((Particles[i].GetName() == "Pi-")&&(Particles[j].GetName()=="K-"))||((Particles[i].GetName() == "K-")&&(Particles[j].GetName()=="Pi-"))){
                        HistoInvMassPionKaonConcordant->Fill(inv_mass);
                        }
                    } else if ((Particles[i].GetCharge()*Particles[j].GetCharge() < 0)){
                        HistoInvMassDiscordi->Fill(inv_mass);
                        if(((Particles[i].GetName() == "Pi+")&&(Particles[j].GetName()=="K-"))||((Particles[i].GetName() == "K-")&&(Particles[j].GetName()=="Pi+"))){
                        HistoInvMassPionKaonDiscordant->Fill(inv_mass);
                        }
                        if(((Particles[i].GetName() == "Pi-")&&(Particles[j].GetName()=="K+"))||((Particles[i].GetName() == "K+")&&(Particles[j].GetName()=="Pi-"))){
                        HistoInvMassPionKaonDiscordant->Fill(inv_mass);
                        }
                    }
                }
             }
        }

        //printing status of the program
        if (i%1000==0){
            std::cout << "Progressivo: " << i/1000 << "%\n";
        }


    }

    //drawing the histogram on the canva
    Canvas->cd(1);
    HistoInvMassPionKaonDiscordant->Draw();
    Canvas->cd(2);
    HistoInvMassPionKaonConcordant->Draw();
    Canvas->cd(3);
    HistoInvMassDecayed->Draw();
    Canvas->cd(4);
    HistoTransverseMomentum->Draw();
    Canvas->cd(5);
    HistoEnergy->Draw();
    Canvas->cd(6);
    HistoInvMass->Draw();
    Canvas->cd(7);
    HistoInvMassConcordi->Draw();
    Canvas->cd(8);
    HistoInvMassDiscordi->Draw();


    File->Write();
    File->Close();
    gBenchmark->Show("");
}

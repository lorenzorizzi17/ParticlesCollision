#include"Particle.hpp"
#include"TRandom.h"
#include"TH1F.h"
#include<vector>
#include"TCanvas.h"

void Draw() {
    
    Particle::AddParticleType("Pi+",0.13957,1,0);
    Particle::AddParticleType("Pi-",0.13957,-1,0);
    Particle::AddParticleType("K+", 0.49367,1,0);
    Particle::AddParticleType("K-", 0.49367,-1,0);
    Particle::AddParticleType("P+",0.93827,1,0);
    Particle::AddParticleType("P-",0.93827,-1,0);
    Particle::AddParticleType("K*",0.89166,0,0.050);

    gRandom->SetSeed();
    TCanvas* Canvas = new TCanvas("c","c", 1400,800);
    Canvas->Divide(4,2);
    

    TH1F* HistoPhi = new TH1F("Azimuthal angle", "Azimuthal angle", 10000,0,2*M_PI);
    TH1F* HistoTheta = new TH1F("Polar angle", "Polar angle", 10000,0,M_PI);
    TH1F* HistoParticleType = new TH1F("Particle type", "Particle type",7,0,7);
    TH1F* HistoMomentum = new TH1F("Momentum distribution", "Momentum distribution", 1000,0,6);
    TH1F* HistoTransverseMomentum = new TH1F("Transverse momentum", "Transverse momentum",1000,0,4);
    TH1F* HistoEnergy = new TH1F("Energy distribution", "Energy distribution", 1000, 0, 3);
    TH1F* HistoInvMass = new TH1F("Invariant mass distribution", "Invariant mass distribution", 10000,0,6);
    TH1F* HistoInvMassConcordi = new TH1F("Inv. mass distribution concordi", "Inv. mass distribution concordi",10000,0,6);
    TH1F* HistoInvMassDiscordi = new TH1F("Inv. mass distribution discordi", "Inv. mass distribution discordi",10000,0,6);
    TH1F* HistoInvMassPionPKaonN = new TH1F("Inv. mass distr. Pi+/K-", "Inv. mass distr. Pi+/K-", 10000,0,6);
    TH1F* HistoInvMassPionNKaonP = new TH1F("Inv. mass distr. Pi-/K+", "Inv. mass distr. Pi-/K+", 10000,0,6);
    TH1F* HistoInvMassDecayed = new TH1F("Inv. mass distr. decayed", "Inv. mass distr. decayed",1000,0.5,1.5);
    std::array<Particle,120> Particles;
    for (int i = 0; i < 1E4; i++)
    {
        Particle* LastInstance;
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
                HistoPhi->Fill(phi);
                HistoMomentum->Fill(pmod);
                HistoTransverseMomentum->Fill(std::sqrt(px*px+py*py));
                HistoEnergy->Fill(particle.GetEnergy());
            } else if (test<0.8) {
                particle.SetIndex("Pi-");
                HistoParticleType->Fill(1);
                HistoPhi->Fill(phi);
                HistoMomentum->Fill(pmod);
                HistoTransverseMomentum->Fill(std::sqrt(px*px+py*py));
                HistoEnergy->Fill(particle.GetEnergy());
            } else if (test < 0.85) {
                particle.SetIndex("K+");
                HistoParticleType->Fill(2);
                HistoPhi->Fill(phi);
                HistoMomentum->Fill(pmod);
                HistoTransverseMomentum->Fill(std::sqrt(px*px+py*py));
                HistoEnergy->Fill(particle.GetEnergy());
            } else if (test < 0.9) {
                particle.SetIndex("K-");
                HistoParticleType->Fill(3);
                HistoPhi->Fill(phi);
                HistoMomentum->Fill(pmod);
                HistoTransverseMomentum->Fill(std::sqrt(px*px+py*py));
                HistoEnergy->Fill(particle.GetEnergy());
            } else if (test <0.945) {
                particle.SetIndex("P+");
                HistoParticleType->Fill(4);
                HistoPhi->Fill(phi);
                HistoMomentum->Fill(pmod);
                HistoTransverseMomentum->Fill(std::sqrt(px*px+py*py));
                HistoEnergy->Fill(particle.GetEnergy());
            } else if (test < 0.99) {
                particle.SetIndex("P-");
                HistoParticleType->Fill(5);
                HistoPhi->Fill(phi);
                HistoMomentum->Fill(pmod);
                HistoTransverseMomentum->Fill(std::sqrt(px*px+py*py));
                HistoEnergy->Fill(particle.GetEnergy());
            } else {
                particle.SetIndex("K*");
                HistoParticleType->Fill(6);
            }
            Particles[j] = particle;
            LastInstance = &particle;
        }
        int k = 100;
        //checking for K* particles
        for (int i{0}; i < 100; i++) {
            if (Particles[i].GetName() == "K*") {
                //making the K*  particle decay
                double test2 = gRandom->Rndm();
                if (test2 < 0.5) {
                    Particle PioneP = Particle("Pi+");
                    Particle KaoneN = Particle("K-");
                    Particles[i].Decay2body(PioneP, KaoneN);
                    Particles[k] = PioneP;
                    Particles[k+1] = KaoneN; 
                    k = k+2;
                    double PionePMomentum = sqrt(PioneP.GetPx()*PioneP.GetPx()+PioneP.GetPy()*PioneP.GetPy()+PioneP.GetPz()*PioneP.GetPz());
                    double KaoneNMomentum = sqrt(KaoneN.GetPx()*KaoneN.GetPx()+KaoneN.GetPy()*KaoneN.GetPy()+KaoneN.GetPz()*KaoneN.GetPz());
                    double PionePTransverseMomentum = sqrt(PioneP.GetPx()*PioneP.GetPx()+PioneP.GetPy()*PioneP.GetPy());
                    double KaoneNTransverseMomentum = sqrt(KaoneN.GetPx()*KaoneN.GetPx()+KaoneN.GetPy()*KaoneN.GetPy());
                    HistoMomentum->Fill(PionePMomentum);
                    HistoMomentum->Fill(KaoneNMomentum);
                    HistoTransverseMomentum->Fill(PionePTransverseMomentum);
                    HistoTransverseMomentum->Fill(KaoneNTransverseMomentum);
                    HistoEnergy->Fill(PioneP.GetEnergy());
                    HistoEnergy->Fill(KaoneN.GetEnergy());
                    HistoInvMassDecayed->Fill(PioneP.InvMass(KaoneN));
                } else {
                    Particle PioneN = Particle("Pi-");
                    Particle KaoneP = Particle("K+");
                    Particles[i].Decay2body(PioneN, KaoneP);
                    Particles[k] = PioneN;
                    Particles[k+1] = KaoneP; 
                    k = k+2;
                    double PioneNMomentum = sqrt(PioneN.GetPx()*PioneN.GetPx()+PioneN.GetPy()*PioneN.GetPy()+PioneN.GetPz()*PioneN.GetPz());
                    double KaonePMomentum = sqrt(KaoneP.GetPx()*KaoneP.GetPx()+KaoneP.GetPy()*KaoneP.GetPy()+KaoneP.GetPz()*KaoneP.GetPz());
                    double PioneNTransverseMomentum = sqrt(PioneN.GetPx()*PioneN.GetPx()+PioneN.GetPy()*PioneN.GetPy());
                    double KaonePTransverseMomentum = sqrt(KaoneP.GetPx()*KaoneP.GetPx()+KaoneP.GetPy()*KaoneP.GetPy());
                    HistoMomentum->Fill(PioneNMomentum);
                    HistoMomentum->Fill(KaonePMomentum);
                    HistoTransverseMomentum->Fill(PioneNTransverseMomentum);
                    HistoTransverseMomentum->Fill(KaonePTransverseMomentum);
                    HistoEnergy->Fill(PioneN.GetEnergy());
                    HistoEnergy->Fill(KaoneP.GetEnergy());
                    HistoInvMassDecayed->Fill(PioneN.InvMass(KaoneP));
                }
            }
        }
        //computing the inv mass
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < k; j++)
            {
                if ((i!=j)&&(Particles[i].GetName() != "K*")&&(Particles[j].GetName() != "K*")) {
                    double inv_mass = Particles[i].InvMass(Particles[j]);
                    HistoInvMass->Fill(inv_mass);
                    if (Particles[i].GetCharge() *Particles[j].GetCharge() > 0){
                        HistoInvMassConcordi->Fill(inv_mass);
                    } else {
                        HistoInvMassDiscordi->Fill(inv_mass);
                    }

                    if((Particles[i].GetName() == "Pi+")&&(Particles[j].GetName()=="K-")){
                        HistoInvMassPionPKaonN->Fill(inv_mass);
                    }
                    if((Particles[i].GetName() == "Pi-")&&(Particles[j].GetName()=="K+")){
                        HistoInvMassPionNKaonP->Fill(inv_mass);
                    }
                }
            }
        }
        


    }
    Canvas->cd(1);
    HistoInvMassPionPKaonN->Draw();
    Canvas->cd(2);
    HistoInvMassPionNKaonP->Draw();
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
}
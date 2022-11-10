#include"Particle.hpp"
#include"TRandom.h"
#include"TH1F.h"
#include<vector>
#include"TCanvas.h"
#include"TFile.h"
#include"TF1.h"




void DoAnalysis() {
    TFile* Histograms = new TFile("Particles.root");
    TH1F* HistoParticleType = (TH1F*)Histograms->Get("Particle type");
    TH1F* HistoPolarAngle = (TH1F*)Histograms->Get("Polar angle");
    TH1F* HistoAzimuthalAngle = (TH1F*)Histograms->Get("Azimuthal angle");
    TH1F* HistoMomentum = (TH1F*)Histograms->Get("Momentum distribution");
    TH1F* HistoTransverseMomentum = (TH1F*)Histograms->Get("Transverse momentum");
    TH1F* HistoEnergy = (TH1F*)Histograms->Get("Energy distribution");
    TH1F* HistoInvMass = (TH1F*)Histograms->Get("Invariant mass distribution");
    TH1F* HistoInvMassConcordi = (TH1F*)Histograms->Get("Inv. mass distribution concordi");
    TH1F* HistoInvMassDiscordi = (TH1F*)Histograms->Get("Inv. mass distribution discordi");
    TH1F* HistoInvMassPionPKaonN = (TH1F*)Histograms->Get("Inv. mass distr. Pi+/K-");
    TH1F* HistoInvMassPionNKaonP = (TH1F*)Histograms->Get("Inv. mass distr. Pi-/K+");
    TH1F* HistoInvMassDecayed = (TH1F*)Histograms->Get("Inv. mass distr. decayed");
    std::cout << "*************************************" << '\n' <<
       "NUMBERS OF ENTRIES IN HISTOGRAMS: \n" <<
       "Particle Type Histogram: " << HistoParticleType->GetEntries() << '\n' <<
       "Polar Angle Histogram: " << HistoPolarAngle->GetEntries() << '\n' <<
       "Azimuthal Angle Histogram: " << HistoAzimuthalAngle->GetEntries() << '\n'<<
       "Total momentum Histogram: " << HistoMomentum->GetEntries() << '\n'<<
       "Transverse momentum Histogram: " << HistoTransverseMomentum->GetEntries() << '\n'<<
       "Energy Histogram: " << HistoEnergy->GetEntries() << '\n'<<
       "Total invariant Mass Histogram: " << HistoInvMass->GetEntries() << '\n'<<
       "Invariant mass of concordant particles histogram: " << HistoInvMassConcordi->GetEntries() << '\n' <<
       "Invariant mass of discordant particles histogram: " << HistoInvMassDiscordi->GetEntries() << '\n' <<
       //"Invariant mass of Pi+ and K- particles histogram: " << HistoInvMassPionPKaonN->GetEntries() << '\n' <<
       //"Invariant mass of Pi- and K+ particles histogram: " << HistoInvMassPionNKaonP->GetEntries() << '\n' <<
       "Invariant mass of decayed particles histogram: " << HistoInvMassDecayed->GetEntries() << '\n' <<
        "*************************************\n" ;
    std::cout << "*************************************" << '\n' <<
    "PROPORTIONS OF GENERATED PARTICLES: '\n'" <<
    "N. of protons generated: " << HistoParticleType->GetBinContent(1) << " +/- " << HistoParticleType->GetBinError(1) << '\n'<<
    "N. of antiprotons generated: " << HistoParticleType->GetBinContent(2) << " +/- " << HistoParticleType->GetBinError(2) << '\n' <<
    "N. of pions generated: " << HistoParticleType->GetBinContent(3) << " +/- " << HistoParticleType->GetBinError(3) << '\n' <<
    "N. of antipions generated: " << HistoParticleType->GetBinContent(4) << " +/- " << HistoParticleType->GetBinError(4) << '\n' <<
    "N. of kaons generated: " << HistoParticleType->GetBinContent(5) << " +/- " << HistoParticleType->GetBinError(5) << '\n' <<
    "N. of antikaons generated: " << HistoParticleType->GetBinContent(6) << " +/- " << HistoParticleType->GetBinError(6) << '\n' <<
    "N. of kaons * generated: " << HistoParticleType->GetBinContent(7) << " +/- " << HistoParticleType->GetBinError(7) << '\n' <<
    "*************************************\n" ;


    TF1* UniformFunctionPolar = new TF1("Uniform Function Polar", "[0]", 0, M_PI);
    TF1* UniformFunctionAzimuthal = new TF1("Uniform Function Azimuthal", "[0]",0,2*M_PI);
    HistoAzimuthalAngle->Fit(UniformFunctionAzimuthal, "QQQQ");
    HistoPolarAngle->Fit(UniformFunctionPolar, "QQQQ");
    std::cout << "*************************************" << '\n' <<
    "Azimuthal angle fitting: " << UniformFunctionAzimuthal->GetParameter(0) << " +/- " << UniformFunctionAzimuthal->GetParError(0) <<
     "                chi2/NDF = " << UniformFunctionAzimuthal->GetChisquare()/UniformFunctionAzimuthal->GetNDF() << '\n' <<
    "Polar angle fitting: " << UniformFunctionPolar->GetParameter(0) << " +/- " << UniformFunctionPolar->GetParError(0) <<
    "                chi2/NDF = " << UniformFunctionPolar->GetChisquare()/UniformFunctionPolar->GetNDF() << '\n' <<
    "*************************************" << '\n' ;


    TF1* ExpoFunction = new TF1("Expo function", "[0]*exp(-[1]*x)",0,3);
    HistoMomentum->Fit(ExpoFunction, "Q");
    std::cout << "*************************************" << '\n' <<
    "Mean momentum of generated particles: " << ExpoFunction->GetParameter(1) << '\n' <<
    "*************************************" << '\n';

}
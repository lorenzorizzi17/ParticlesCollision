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
    TH1F* HistoTransverseMomentum = (TH1F*)Histograms->Get("Transverse momentum distribution");
    TH1F* HistoEnergy = (TH1F*)Histograms->Get("Energy distribution");
    TH1F* HistoInvMass = (TH1F*)Histograms->Get("Invariant mass distribution");
    TH1F* HistoInvMassConcordant = (TH1F*)Histograms->Get("Inv. mass distr. concordant particles");
    //HistoInvMassConcordant->Sumw2();
    TH1F* HistoInvMassDiscordant = (TH1F*)Histograms->Get("Inv. mass distr. discordant");
    //HistoInvMassDiscordant->Sumw2();
    TH1F* HistoInvMassPionPKaonN = (TH1F*)Histograms->Get("Inv. mass distr. concordant Pi and K particles");
    TH1F* HistoInvMassPionNKaonP = (TH1F*)Histograms->Get("Inv. mass distr. discordant Pi and K particles");
    TH1F* HistoInvMassDecayed = (TH1F*)Histograms->Get("Inv. mass distr. daughter particles");
    std::cout << "*************************************" << '\n' <<
       "NUMBERS OF ENTRIES IN HISTOGRAMS: \n" <<
       "Particle Type Histogram: " << HistoParticleType->GetEntries() << '\n' <<
       "Polar Angle Histogram: " << HistoPolarAngle->GetEntries() << '\n' <<
       "Azimuthal Angle Histogram: " << HistoAzimuthalAngle->GetEntries() << '\n'<<
       "Total momentum Histogram: " << HistoMomentum->GetEntries() << '\n'<<
       "Transverse momentum Histogram: " << HistoTransverseMomentum->GetEntries() << '\n'<<
       "Energy Histogram: " << HistoEnergy->GetEntries() << '\n'<<
       "Total invariant Mass Histogram: " << HistoInvMass->GetEntries() << '\n'<<
       "Invariant mass of concordant particles histogram: " << HistoInvMassConcordant->GetEntries() << '\n' <<
       "Invariant mass of discordant particles histogram: " << HistoInvMassDiscordant->GetEntries() << '\n' <<
       "Invariant mass of Pi+ and K- particles histogram: " << HistoInvMassPionPKaonN->GetEntries() << '\n' <<
       "Invariant mass of Pi- and K+ particles histogram: " << HistoInvMassPionNKaonP->GetEntries() << '\n' <<
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
    "Mean momentum of generated particles: " << ExpoFunction->GetParameter(1) <<
     "                chi2/NDF = " << ExpoFunction->GetChisquare()/ExpoFunction->GetNDF() << '\n' <<
    "*************************************" << '\n';


    TF1* Gaus = new TF1("Gaus", "[0]*exp(-((x-[1])^2)/(2*[2]^2))",0,3);
    Gaus->SetParameter(0,1000);
    Gaus->SetParameter(1,0.8);
    Gaus->SetParameter(2,0.5);
    TCanvas* Canvas = new TCanvas("Canvas", "Canvas", 800,600);
    Canvas->Divide(3,2);
    Canvas->cd(1);
    HistoInvMassConcordant->GetYaxis()->SetTitle("Invariant mass (GeV/c^2)");
    HistoInvMassConcordant->GetXaxis()->SetTitle("Entries");
    HistoInvMassConcordant->SetFillColor(kBlue);
    HistoInvMassConcordant->DrawCopy("HIST");
    Canvas->cd(2);
    HistoInvMassDiscordant->GetYaxis()->SetTitle("Invariant mass (GeV/c^2)");
    HistoInvMassDiscordant->GetXaxis()->SetTitle("Entries");
    HistoInvMassDiscordant->SetFillColor(kBlue);
    HistoInvMassDiscordant->DrawCopy("HIST");
    Canvas->cd(3);
    TH1F* HistoDifference = new TH1F("Difference histogram, concordant/discordant", "Difference histogram, discordant-concordant",160,0,3);
    HistoDifference->Add(HistoInvMassDiscordant,HistoInvMassConcordant,1,-1);
    HistoDifference->Fit(Gaus, "Q");
    std::cout << "*************************************" << '\n' <<
    "K* mass according to first histograms difference: (" << Gaus->GetParameter(1) << " +/- " << Gaus->GetParError(1) << ") GeV/c^2 \n" <<
    "K* lenght according to first histograms difference: (" << Gaus->GetParameter(2) << " +/- " << Gaus->GetParError(2) << ") GeV/c^2 \n" <<
     "                chi2/NDF = " << Gaus->GetChisquare()/Gaus->GetNDF() << '\n' <<
    "*************************************" << '\n';
    HistoDifference->GetYaxis()->SetTitle("Invariant mass (GeV/c^2)");
    HistoDifference->GetXaxis()->SetTitle("Entries");
    HistoDifference->SetFillColor(kBlue);
    HistoDifference->DrawCopy("HIST");
    HistoDifference->DrawCopy("SAME,FUNC");
    Canvas->cd(4);
    HistoInvMassPionPKaonN->GetYaxis()->SetTitle("Invariant mass (GeV/c^2)");
    HistoInvMassPionPKaonN->GetXaxis()->SetTitle("Entries");
    HistoInvMassPionPKaonN->SetFillColor(kBlue);
    HistoInvMassPionPKaonN->DrawCopy("HIST");
    Canvas->cd(5);
    HistoInvMassPionNKaonP->GetYaxis()->SetTitle("Invariant mass (GeV/c^2)");
    HistoInvMassPionNKaonP->GetXaxis()->SetTitle("Entries");
    HistoInvMassPionNKaonP->SetFillColor(kBlue);
    HistoInvMassPionNKaonP->DrawCopy("HIST");
    Canvas->cd(6);
    TH1F* HistoDifference1 = new TH1F("Difference histogram Pi/K, concordant/discordant", "Difference histogram Pi/K, discordant-concordant",160,0,3);
    HistoDifference1->Add(HistoInvMassPionPKaonN,HistoInvMassPionNKaonP,-1,1);
    HistoDifference1->Fit(Gaus, "Q");
    std::cout << "*************************************" << '\n' <<
    "K* mass according to second histograms difference: (" << Gaus->GetParameter(1) << " +/- " << Gaus->GetParError(1) << ") GeV/c^2 \n" <<
    "K* lenght according to second histograms difference: (" << Gaus->GetParameter(2) << " +/- " << Gaus->GetParError(2) << ") GeV/c^2 \n" <<
     "                chi2/NDF = " << Gaus->GetChisquare()/Gaus->GetNDF() << '\n' <<
    "*************************************" << '\n';
    HistoDifference1->GetXaxis()->SetTitle("Invariant mass (GeV/c^2)");
    HistoDifference1->GetXaxis()->SetTitle("Entries");
    HistoDifference1->SetFillColor(kBlue);
    HistoDifference1->DrawCopy("HIST");
    HistoDifference1->DrawCopy("SAME,FUNC");
}

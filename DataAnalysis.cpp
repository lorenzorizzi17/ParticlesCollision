#include"Particle.hpp"
#include"TRandom.h"
#include"TH1F.h"
#include<vector>
#include"TCanvas.h"
#include"TFile.h"
#include"TF1.h"

void setStyle(TH1F* h, const char* title, const char* XaxisName){
    h->SetTitle(title);
    h->GetYaxis()->SetTitle("Entries");
    h->GetXaxis()->SetTitle(XaxisName);
    h->SetFillColor(33);
    h->SetLineColor(4);
}

void analysis() {
    TFile* Histograms = new TFile("Particles.root");
    //getting the histograms from root file 
    TH1F* HistoParticleType = (TH1F*)Histograms->Get("Particle type");
    TH1F* HistoPolarAngle = (TH1F*)Histograms->Get("Polar angle");
    TH1F* HistoAzimuthalAngle = (TH1F*)Histograms->Get("Azimuthal angle");
    TH1F* HistoMomentum = (TH1F*)Histograms->Get("Momentum distribution");
    TH1F* HistoTransverseMomentum = (TH1F*)Histograms->Get("Transverse momentum distribution");
    TH1F* HistoEnergy = (TH1F*)Histograms->Get("Energy distribution");
    TH1F* HistoInvMass = (TH1F*)Histograms->Get("Invariant mass distribution");
    TH1F* HistoInvMassConcordant = (TH1F*)Histograms->Get("Inv. mass distr. concordant particles");
    TH1F* HistoInvMassDiscordant = (TH1F*)Histograms->Get("Inv. mass distr. discordant particles");
    TH1F* HistoInvMassPionPKaonN = (TH1F*)Histograms->Get("Inv. mass distr. concordant Pi and K particles");
    TH1F* HistoInvMassPionNKaonP = (TH1F*)Histograms->Get("Inv. mass distr. discordant Pi and K particles");
    TH1F* HistoInvMassDecayed = (TH1F*)Histograms->Get("Inv. mass distr. daughter particles");

    //printing entries
    std::cout << "**************************************************************************" << '\n' <<
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
        "**************************************************************************\n" ;
    std::cout << "**************************************************************************" << '\n' <<
    "PROPORTIONS OF GENERATED PARTICLES: '\n'" <<
    "N. of pions generated: " << HistoParticleType->GetBinContent(1) << " +/- " << HistoParticleType->GetBinError(1) << '\n'<<
    "N. of antipions generated: " << HistoParticleType->GetBinContent(2) << " +/- " << HistoParticleType->GetBinError(2) << '\n' <<
    "N. of protons generated: " << HistoParticleType->GetBinContent(3) << " +/- " << HistoParticleType->GetBinError(3) << '\n' <<
    "N. of antiprotons generated: " << HistoParticleType->GetBinContent(4) << " +/- " << HistoParticleType->GetBinError(4) << '\n' <<
    "N. of kaons generated: " << HistoParticleType->GetBinContent(5) << " +/- " << HistoParticleType->GetBinError(5) << '\n' <<
    "N. of antikaons generated: " << HistoParticleType->GetBinContent(6) << " +/- " << HistoParticleType->GetBinError(6) << '\n' <<
    "N. of kaons * generated: " << HistoParticleType->GetBinContent(7) << " +/- " << HistoParticleType->GetBinError(7) << '\n' <<
    "**************************************************************************\n" ;


    TF1* UniformFunctionPolar = new TF1("Uniform Function Polar", "[0]", 0, M_PI);
    TF1* UniformFunctionAzimuthal = new TF1("Uniform Function Azimuthal", "[0]",0,2*M_PI);
    HistoAzimuthalAngle->Fit(UniformFunctionAzimuthal, "Q");
    HistoPolarAngle->Fit(UniformFunctionPolar, "Q");
    std::cout << "**************************************************************************" << '\n' <<
    "Azimuthal angle fitting: " << UniformFunctionAzimuthal->GetParameter(0) << " +/- " << UniformFunctionAzimuthal->GetParError(0) <<
     "                chi2/NDF = " << UniformFunctionAzimuthal->GetChisquare()/UniformFunctionAzimuthal->GetNDF()  << ", probability: " << UniformFunctionAzimuthal->GetProb() << '\n' <<
    "Polar angle fitting: " << UniformFunctionPolar->GetParameter(0) << " +/- " << UniformFunctionPolar->GetParError(0) <<
    "                chi2/NDF = " << UniformFunctionPolar->GetChisquare()/UniformFunctionPolar->GetNDF() << ", probability: " << UniformFunctionPolar->GetProb() << '\n' <<
    "**************************************************************************" << '\n' ;


    TF1* ExpoFunction = new TF1("Expo function", "[0]*exp(-[1]*x)",0,3);
    HistoMomentum->Fit(ExpoFunction, "QQQ");
    std::cout << "**************************************************************************" << '\n' <<
    "Mean momentum of generated particles: " << ExpoFunction->GetParameter(1) <<
     "                chi2/NDF = " << ExpoFunction->GetChisquare()/ExpoFunction->GetNDF() << ", probability: " << ExpoFunction->GetProb() << '\n' <<
    "**************************************************************************" << '\n';


    TF1* Gaus = new TF1("Gaus", "[0]*exp(-((x-[1])^2)/(2*[2]^2))",0.5,1.1);
    Gaus->SetParameter(0,17000);
    Gaus->SetParameter(1,0.9);
    Gaus->SetParameter(2,0.5);

    TCanvas* c1 = new TCanvas("Invariant mass of distribution concordant particles","Invariant mass distribution of concordant particles",800,600);
    setStyle(HistoInvMassConcordant,"Invariant mass distribution of concordant particles","Invariant mass (GeV/c^2)");
    c1->Print("Invariant mass of concordant particles.pdf");
    c1->Print("Invariant mass of concordant particles.root");
    c1->Print("Invariant mass of concordant particles.C");

    TCanvas* c2 = new TCanvas("Invariant mass distribution of discordant particles", "Invariant mass of discordant particles",800,600);
    setStyle(HistoInvMassDiscordant,"Invariant mass distribution of discordant particles","Invariant mass (GeV/c^2)");
    HistoInvMassDiscordant->Draw("HIST");
    c2->Print("Invariant mass of discordant particles.pdf");
    c2->Print("Invariant mass of discordant particles.root");
    c2->Print("Invariant mass of discordant particles.C");

    TCanvas* c3 = new TCanvas("Difference histogram, concordant-discordant", "Difference histogram, concordant-discordant", 800,600);
    TH1F* HistoDifference = new TH1F("Difference histogram, concordant/discordant", "Difference histogram, discordant-concordant",150,0,3);
    HistoDifference->Add(HistoInvMassDiscordant,HistoInvMassConcordant,1,-1);
    HistoDifference->Fit(Gaus, "Q");
    std::cout << "**************************************************************************" << '\n' <<
    "K* mass according to first histograms difference: (" << Gaus->GetParameter(1) << " +/- " << Gaus->GetParError(1) << ") GeV/c^2 \n" <<
    "K* lenght according to first histograms difference: (" << Gaus->GetParameter(2) << " +/- " << Gaus->GetParError(2) << ") GeV/c^2 \n" <<
     "                chi2/NDF = " << Gaus->GetChisquare()/Gaus->GetNDF() <<  ", probability: " << Gaus->GetProb() << '\n' <<
    "**************************************************************************" << '\n';
    setStyle(HistoDifference,"Difference histogram, concordant-discordant particles","Invariant mass (GeV/c^2)");
    HistoDifference->Draw("HIST");
    HistoDifference->Draw("SAME,FUNC");
    c3->Print("Difference histogram, concordant-discordant.pdf");
    c3->Print("Difference histogram, concordant-discordant.root");
    c3->Print("Difference histogram, concordant-discordant.C");


    TCanvas* c4 = new TCanvas("Invariant mass of #pi-k concordant particles", "Invariant mass of #pi-k concordant particles",800,600);
    setStyle(HistoInvMassPionPKaonN, "Invariant mass distribution of concordant #pi-K particles", "Invariant mass (GeV/c^2)" );
    HistoInvMassPionPKaonN->Draw("HIST");
    c4->Print("Invariant mass of #pi-k concordant particles.pdf");
    c4->Print("Invariant mass of #pi-k concordant particles.root");
    c4->Print("Invariant mass of #pi-k concordant particles.C");

    TCanvas* c5 = new TCanvas("Invariant mass of #pi/k discordant particles", "Invariant mass of #pi-k discordant particles",800,600);
    setStyle(HistoInvMassPionNKaonP, "Invariant mass distribution of discordant #pi-K particles", "Invariant mass (GeV/c^2)" );
    HistoInvMassPionNKaonP->Draw("HIST");
    c5->Print("Invariant mass of #pi-k discordant particle.pdf");
    c5->Print("Invariant mass of #pi-k discordant particle.root");
    c5->Print("Invariant mass of #pi-k discordant particle.C");

    TCanvas* c6 = new TCanvas("Difference histogram Pi-K, concordant-discordant","Difference histogram Pi-K, concordant-discordant",800,600);
    TH1F* HistoDifference1 = new TH1F("Difference histogram Pi-K, concordant-discordant", "Difference histogram Pi-K, discordant-concordant",150,0,3);
    HistoDifference1->Add(HistoInvMassPionPKaonN,HistoInvMassPionNKaonP,-1,1);
    HistoDifference1->Fit(Gaus, "QR");
    std::cout << "**************************************************************************" << '\n' <<
    "K* mass according to second histograms difference: (" << Gaus->GetParameter(1) << " +/- " << Gaus->GetParError(1) << ") GeV/c^2 \n" <<
    "K* lenght according to second histograms difference: (" << Gaus->GetParameter(2) << " +/- " << Gaus->GetParError(2) << ") GeV/c^2 \n" <<
     "                chi2/NDF = " << Gaus->GetChisquare()/Gaus->GetNDF() << ", probability: " << Gaus->GetProb() << '\n' <<
    "**************************************************************************" << '\n';
    setStyle(HistoDifference1, "Difference histogram, concordant-discordant #pi-K particles","Invariant mass (GeV/c^2)" );
    HistoDifference1->Draw("HIST");
    HistoDifference1->Draw("SAME,FUNC");
    c6->Print("Difference histogram Pi-K, concordant-discordant.pdf");
    c6->Print("Difference histogram Pi-K, concordant-discordant.C");
    c6->Print("Difference histogram Pi-K, concordant-discordant.root");

    TCanvas* c7 = new TCanvas("Invariant mass of daughter particles coupled","Invariant mass of daughter particles coupled",800,600);
    TF1* Gaus2 = new TF1("Gaus2", "[0]*exp(-((x-[1])^2)/(2*[2]^2))",0.5,1.1);
    Gaus2->SetParameter(0,12000);
    Gaus2->SetParameter(1,0.8);
    Gaus2->SetParameter(2,0.5);
    HistoInvMassDecayed->GetXaxis()->SetTitle("Invariant mass (GeV/c^2)");
    HistoInvMassDecayed->GetXaxis()->SetTitle("Entries");
    HistoInvMassDecayed->SetFillColor(kBlue);
    HistoInvMassDecayed->Fit(Gaus2, "Q");
    std::cout << "**************************************************************************" << '\n' <<
    "K* mass according to daughter particles histograms: (" << Gaus2->GetParameter(1) << " +/- " << Gaus2->GetParError(1) << ") GeV/c^2 \n" <<
    "K* lenght according to daughter particles histograms: (" << Gaus2->GetParameter(2) << " +/- " << Gaus2->GetParError(2) << ") GeV/c^2 \n" <<
     "                chi2/NDF = " << Gaus2->GetChisquare()/Gaus2->GetNDF() << ", probability: " << Gaus2->GetProb() << '\n' <<
    "**************************************************************************" << '\n';
    setStyle(HistoInvMassDecayed, "Invariant mass distribution of daughter particles", "Invariant mass (GeV/c^2)");
    HistoInvMassDecayed->Draw("HIST");
    HistoInvMassDecayed->Draw("SAME,FUNC");
    c7->Print("Difference histogram Pi-K, concordant-discordant.pdf");
    c7->Print("Difference histogram Pi-K, concordant-discordant.C");
    c7->Print("Difference histogram Pi-K, concordant-discordant.root");
}

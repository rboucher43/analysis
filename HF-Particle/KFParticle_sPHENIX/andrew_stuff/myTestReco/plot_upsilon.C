#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

void plot_upsilon(const char* filename="/sphenix/user/rboucher43/analysis/HF-Particle/KFParticle_sPHENIX/andrew_stuff/myTestReco/10000_events_upsilon_andrew.root") 
    {
    // Open the ROOT file
    TFile *file = new TFile(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }
    
    // Access the decay tree
    TTree *decayTree = (TTree*)file->Get("DecayTree");
    if (!decayTree) 
      {
        std::cerr << "Error: Unable to find DecayTree in file " << filename << std::endl;
        file->Close();
        return;
    }
    
    // Create histogram
    TH1F *hist = new TH1F("InvariantMassHist", "Invariant Mass Distribution", 100, 0, 10);
    
    // Fill histogram with data from the tree
    Float_t invariantMass;
    decayTree->SetBranchAddress("Upsilon_mass", &invariantMass);
    Long64_t nEntries = decayTree->GetEntries();
    std::cout << "nEntries "<<nEntries << std::endl;
    for (Long64_t i = 0; i < nEntries; ++i) {
        decayTree->GetEntry(i);
	std::cout << "invariant_mass "<< invariantMass << std::endl;
        hist->Fill(invariantMass);
    }
    
    // Create canvas and draw histogram
    TCanvas *canvas = new TCanvas("canvas", "Upsilon Mass", 800, 600);
    hist->DrawCopy();
    //canvas->Draw();
    
    // Cleanup
    file->Close();
}

#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"

// from ND280 software
#include "TNRooTrackerVtx.hxx"

#define DEBUG(x) std::cout << x << std::endl

class NRooTrackerConverter {

 public:
  
  NRooTrackerConverter(std::string inputFileName, std::string outputFileName);
  void SetAddresses();
  void FillVtx(int entry);
  inline int GetEntries() { return fRooTrackerTree->GetEntries(); };
  inline void Fill() { fOutputTree->Fill(); }; 
  inline void Write() { fOutputTree->Write(); };
  inline void Close() { fOutputFile->Close(); };

 private:

  ND::NRooTrackerVtx *fCurrNeutVtx;
  TTree *fRooTrackerTree;
  TTree *fOutputTree;
  TFile *fFile;
  TFile *fOutputFile;

  Int_t fNVtx;        ///< Number of vertices
  TClonesArray *fVtx; ///< Array of vertex objects
  
  int fOrigTreeEntryNumber;   ///< Entry in original rootracker tree.
  int fInputTreeEntryNumber;   ///< The entry number in the current rootracker file.
  int fOrigInputTreeEntries;   ///< The number of entries in the original input tree. Important for working  out POTs.
  double fOrigInputTreePOT;
  std::string fOrigInputFileName;   ///< The original rootracker tree POT.
  std::string fGeneratorName;   ///< The generator name. Currently genie or neut.
  std::string fInputTreeName;   ///< The input rootracker tree name. Currently gRooTracker or nRooTracker.
  double fTimeInSpill;   ///< The time within the spill for vertex.
  int fTruthVertexID;   ///< The ID of the truth vertex created from this generator vertex.

  bool fHaveNFBranches;
};


NRooTrackerConverter::NRooTrackerConverter(std::string inputFileName, std::string outputFileName) {

  fFile = new TFile(inputFileName.c_str());
  if (!fFile) {
    throw std::runtime_error("Bad input file!");
  }

  fRooTrackerTree = fFile->Get<TTree>("nRooTracker");

  fCurrNeutVtx = new ND::NRooTrackerVtx();
  fVtx = new TClonesArray("ND::NRooTrackerVtx", 100);
  
  fOutputFile = new TFile(outputFileName.c_str(), "NEW");
  if (!fOutputFile) {
    throw std::runtime_error("Bad output file!");
  }

  fOutputFile->cd();
  fOutputTree = new TTree("NRooTrackerVtx", "The Neut information pass-through from nRooTracker tree");

  fOutputTree->Branch("NVtx", &fNVtx, "NVtx/I", 32000);
  fOutputTree->Branch("Vtx", "TClonesArray", &fVtx, 32000, 1);

  DEBUG("NRooTrackerConverter created");

}


void NRooTrackerConverter::FillVtx(int entry) {
  
  fNVtx = 0;
  fCurrNeutVtx->Reset();
  fRooTrackerTree->GetEntry(entry);

  // first fill up all the book-keeping info
  fCurrNeutVtx->OrigFileName->SetString(fOrigInputFileName.c_str());
  fCurrNeutVtx->OrigTreeName->SetString(fInputTreeName.c_str());
  fCurrNeutVtx->OrigEvtNum = fOrigTreeEntryNumber;
  fCurrNeutVtx->OrigTreePOT = fOrigInputTreePOT;
  fCurrNeutVtx->OrigTreeEntries = fOrigInputTreeEntries;
  fCurrNeutVtx->TimeInSpill = fTimeInSpill;
  fCurrNeutVtx->TruthVertexID = fTruthVertexID;

  // Need to convert to the nd280 length units
  fCurrNeutVtx->EvtVtx[0] *= 1000.0;
  fCurrNeutVtx->EvtVtx[1] *= 1000.0;
  fCurrNeutVtx->EvtVtx[2] *= 1000.0;

  ///< Get the name of the generator that generated this event.
  fCurrNeutVtx->GeneratorName->SetString(fGeneratorName.c_str());

  ///< Get the geometry path for this event.
  /// @todo: Make get geometry for IWCD / HK ??????
  fCurrNeutVtx->GeomPath->SetString("BLAAAAA");

  ///< Now create another fRooTracker vtx to add to the vertex
  ///< container.
  ND::NRooTrackerVtx *fNRooTrackerVtx;
  fNRooTrackerVtx = new ((*fVtx)[fNVtx++]) ND::NRooTrackerVtx();

  /// Copy the current neut event. Make sure this is done last so we copy a
  /// complete event.
  fNRooTrackerVtx->Copy(fCurrNeutVtx);

  DEBUG("Vtx filled");
}


void NRooTrackerConverter::SetAddresses() {

  // Now set the tree branch addresses
  fRooTrackerTree->SetBranchAddress("EvtCode", &(fCurrNeutVtx->EvtCode));
  fRooTrackerTree->SetBranchAddress("EvtNum", &(fCurrNeutVtx->EvtNum));
  fRooTrackerTree->SetBranchAddress("EvtXSec", &(fCurrNeutVtx->EvtXSec));
  fRooTrackerTree->SetBranchAddress("EvtDXSec", &(fCurrNeutVtx->EvtDXSec));
  fRooTrackerTree->SetBranchAddress("EvtWght", &(fCurrNeutVtx->EvtWght));
  fRooTrackerTree->SetBranchAddress("EvtProb", &(fCurrNeutVtx->EvtProb));
  fRooTrackerTree->SetBranchAddress("EvtVtx", (fCurrNeutVtx->EvtVtx));
  fRooTrackerTree->SetBranchAddress("StdHepN", &(fCurrNeutVtx->StdHepN));
  fRooTrackerTree->SetBranchAddress("StdHepPdg",
                                    &(fCurrNeutVtx->StdHepPdgTemp));
  fRooTrackerTree->SetBranchAddress("StdHepStatus",
                                    &(fCurrNeutVtx->StdHepStatusTemp));
  fRooTrackerTree->SetBranchAddress("StdHepX4", (fCurrNeutVtx->StdHepX4));
  fRooTrackerTree->SetBranchAddress("StdHepP4", (fCurrNeutVtx->StdHepP4));
  fRooTrackerTree->SetBranchAddress("StdHepPolz", (fCurrNeutVtx->StdHepPolz));
  fRooTrackerTree->SetBranchAddress("StdHepFd", (fCurrNeutVtx->StdHepFdTemp));
  fRooTrackerTree->SetBranchAddress("StdHepLd", (fCurrNeutVtx->StdHepLdTemp));
  fRooTrackerTree->SetBranchAddress("StdHepFm", (fCurrNeutVtx->StdHepFmTemp));
  fRooTrackerTree->SetBranchAddress("StdHepLm", (fCurrNeutVtx->StdHepLmTemp));

  // NEUT > v5.0.7 && MCP > 1 (>10a)
  fRooTrackerTree->SetBranchAddress("NEnvc", &(fCurrNeutVtx->NEnvc));
  fRooTrackerTree->SetBranchAddress("NEipvc", (fCurrNeutVtx->NEipvcTemp));
  fRooTrackerTree->SetBranchAddress("NEpvc", (fCurrNeutVtx->NEpvc));
  fRooTrackerTree->SetBranchAddress("NEiorgvc", (fCurrNeutVtx->NEiorgvcTemp));
  fRooTrackerTree->SetBranchAddress("NEiflgvc", (fCurrNeutVtx->NEiflgvcTemp));
  fRooTrackerTree->SetBranchAddress("NEicrnvc", (fCurrNeutVtx->NEicrnvcTemp));

  fRooTrackerTree->SetBranchAddress("NEnvert", &(fCurrNeutVtx->NEnvert));
  fRooTrackerTree->SetBranchAddress("NEposvert", (fCurrNeutVtx->NEposvert));
  fRooTrackerTree->SetBranchAddress("NEiflgvert",
                                    (fCurrNeutVtx->NEiflgvertTemp));
  fRooTrackerTree->SetBranchAddress("NEnvcvert", &(fCurrNeutVtx->NEnvcvert));
  fRooTrackerTree->SetBranchAddress("NEdirvert", (fCurrNeutVtx->NEdirvert));
  fRooTrackerTree->SetBranchAddress("NEabspvert",
                                    (fCurrNeutVtx->NEabspvertTemp));
  fRooTrackerTree->SetBranchAddress("NEabstpvert",
                                    (fCurrNeutVtx->NEabstpvertTemp));
  fRooTrackerTree->SetBranchAddress("NEipvert", (fCurrNeutVtx->NEipvertTemp));
  fRooTrackerTree->SetBranchAddress("NEiverti", (fCurrNeutVtx->NEivertiTemp));
  fRooTrackerTree->SetBranchAddress("NEivertf", (fCurrNeutVtx->NEivertfTemp));
  // end NEUT > v5.0.7 && MCP > 1 (>10a)

  // NEUT > v5.1.2 && MCP >= 5
  fRooTrackerTree->SetBranchAddress("NEcrsx", &(fCurrNeutVtx->NEcrsx));
  fRooTrackerTree->SetBranchAddress("NEcrsy", &(fCurrNeutVtx->NEcrsy));
  fRooTrackerTree->SetBranchAddress("NEcrsz", &(fCurrNeutVtx->NEcrsz));
  fRooTrackerTree->SetBranchAddress("NEcrsphi", &(fCurrNeutVtx->NEcrsphi));

  // Nucleon FSI information
  fHaveNFBranches = !fRooTrackerTree->SetBranchAddress("NFnvert", &(fCurrNeutVtx->NFnvert));
  if(fHaveNFBranches){
    fRooTrackerTree->SetBranchAddress("NFiflag", (fCurrNeutVtx->NFiflagTEMP));
    fRooTrackerTree->SetBranchAddress("NFx", (fCurrNeutVtx->NFxTEMP));
    fRooTrackerTree->SetBranchAddress("NFy", (fCurrNeutVtx->NFyTEMP));
    fRooTrackerTree->SetBranchAddress("NFz", (fCurrNeutVtx->NFzTEMP));
    fRooTrackerTree->SetBranchAddress("NFpx", (fCurrNeutVtx->NFpxTEMP));
    fRooTrackerTree->SetBranchAddress("NFpy", (fCurrNeutVtx->NFpyTEMP));
    fRooTrackerTree->SetBranchAddress("NFpz", (fCurrNeutVtx->NFpzTEMP));
    fRooTrackerTree->SetBranchAddress("NFe", (fCurrNeutVtx->NFeTEMP));
    fRooTrackerTree->SetBranchAddress("PCascProb", &(fCurrNeutVtx->PCascProb));

    fRooTrackerTree->SetBranchAddress("NFfirststep",
                                      (fCurrNeutVtx->NFfirststepTEMP));
    fRooTrackerTree->SetBranchAddress("NFnstep", &(fCurrNeutVtx->NFnstep));
    fRooTrackerTree->SetBranchAddress("NFecms2", (fCurrNeutVtx->NFecms2TEMP));
    fRooTrackerTree->SetBranchAddress("Prob", (fCurrNeutVtx->ProbTEMP));
    fRooTrackerTree->SetBranchAddress("VertFlagStep", (fCurrNeutVtx->VertFlagStepTEMP));
    fRooTrackerTree->SetBranchAddress("VertFsiRhon", (fCurrNeutVtx->VertFsiRhonTEMP));
    fRooTrackerTree->SetBranchAddress("StepPel", (fCurrNeutVtx->StepPelTEMP));
    fRooTrackerTree->SetBranchAddress("StepPsp", (fCurrNeutVtx->StepPspTEMP));
    fRooTrackerTree->SetBranchAddress("StepPdp", (fCurrNeutVtx->StepPdpTEMP));
  } else {
    std::cout << "Passthrough tree doesn't contain" << std::endl
      << " nucleon FSI tracking information, these branches will remain off." << std::endl;
  }
  fRooTrackerTree->SetBranchAddress("SPIDelta", &(fCurrNeutVtx->SPIDelta));
  fRooTrackerTree->SetBranchAddress("IRadCorrPht", &(fCurrNeutVtx->IRadCorrPht));

  DEBUG("addresses set");
}


int main(int argc, char *argv[]) {
  
  NRooTrackerConverter converter(argv[1], argv[2]);
  converter.SetAddresses();

  for (int entry = 0; entry < converter.GetEntries(); entry++) {
    converter.FillVtx(entry);
    converter.Fill();
  }
  
  converter.Write();
  converter.Close();
}

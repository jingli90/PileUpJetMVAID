//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 17 09:49:57 2015 by ROOT version 6.02/05
// from TTree t/t
// found on file: /eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150709/150709_160648/0000/output_mc_9.root
//////////////////////////////////////////////////////////

#ifndef base_muons_h
#define base_muons_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "Math/GenVector/PtEtaPhiE4D.h"
#include "map"

// Fixed size dimensions of array or collections stored in the TTree if any.
//const Int_t kMaxp4 = 1000;
//const Int_t kMaxgen_p4 = 1000;
#include "base.h" // put the size in to one file to avoid duplication

class base_muons {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   std::vector<float>   *relativeIsoR04_puppiNoMuonWeighted;
   std::vector<float>   *photonIsoR04_puppiNoMuonWeighted;
   std::vector<float>   *neutralHadronIsoR04_puppiNoMuonWeighted;
   std::vector<float>   *relativeIsoR04_puppiWeighted;
   std::vector<float>   *chargedHadronIsoR04_puppiWeighted;
   std::vector<float>   *relativeIsoR04_pfWeighted;
   std::vector<bool>    *isTight;
   std::vector<float>   *neutralHadronIsoR04_puppiWeighted;
   std::vector<bool>    *isHighPt;
   std::vector<bool>    *isSoft;
   std::vector<float>   *photonIsoR03;
   std::vector<float>   *neutralHadronIsoR03;
   //   std::vector<signed char> *gen_charge;
   std::vector<std::map<std::string,bool> > *ids;
   Int_t           p4_;
   Float_t         p4_fCoordinates_fPt[kMaxp4];   //[p4_]
   Float_t         p4_fCoordinates_fEta[kMaxp4];   //[p4_]
   Float_t         p4_fCoordinates_fPhi[kMaxp4];   //[p4_]
   Float_t         p4_fCoordinates_fE[kMaxp4];   //[p4_]
   std::vector<float>   *chargedHadronIsoR04_puppiNoMuonWeighted;
   std::vector<float>   *relativeIsoR03_withEA;
   std::vector<float>   *photonIsoR04_pfWeighted;
   Int_t           gen_p4_;
   Float_t         gen_p4_fCoordinates_fPt[kMaxgen_p4];   //[gen_p4_]
   Float_t         gen_p4_fCoordinates_fEta[kMaxgen_p4];   //[gen_p4_]
   Float_t         gen_p4_fCoordinates_fPhi[kMaxgen_p4];   //[gen_p4_]
   Float_t         gen_p4_fCoordinates_fE[kMaxgen_p4];   //[gen_p4_]
   std::vector<bool>    *has_gen_particle;
   std::vector<float>   *chargedHadronIsoR03;
      //std::vector<signed char> *charge; // test for signed char
	  //std::vector<int8_t> *charge; // test for signed char
	  std::vector<float> *charge; // test for signed char
   std::vector<float>   *photonIsoR04_puppiWeighted;
   std::vector<float>   *y;
   std::vector<bool>    *isLoose;
   std::vector<float>   *neutralHadronIsoR04_pfWeighted;
   std::vector<float>   *gen_y;
   std::vector<float>   *photonIsoR04;
   std::vector<float>   *puChargedHadronIsoR03;
   std::vector<bool>    *isMedium;
   std::vector<float>   *relativeIsoR03_deltaBeta;
   std::vector<float>   *relativeIsoR03;
   std::vector<float>   *chargedHadronIsoR04;
   std::vector<float>   *neutralHadronIsoR04;
   std::vector<float>   *relativeIsoR04_deltaBeta;
   std::vector<float>   *puChargedHadronIsoR04;
   std::vector<float>   *relativeIsoR04;
   std::vector<float>   *relativeIsoR04_withEA;

   // List of branches
   TBranch        *b_relativeIsoR04_puppiNoMuonWeighted;   //!
   TBranch        *b_photonIsoR04_puppiNoMuonWeighted;   //!
   TBranch        *b_neutralHadronIsoR04_puppiNoMuonWeighted;   //!
   TBranch        *b_relativeIsoR04_puppiWeighted;   //!
   TBranch        *b_chargedHadronIsoR04_puppiWeighted;   //!
   TBranch        *b_relativeIsoR04_pfWeighted;   //!
   TBranch        *b_isTight;   //!
   TBranch        *b_neutralHadronIsoR04_puppiWeighted;   //!
   TBranch        *b_isHighPt;   //!
   TBranch        *b_isSoft;   //!
   TBranch        *b_photonIsoR03;   //!
   TBranch        *b_neutralHadronIsoR03;   //!
   //   TBranch        *b_gen_charge;   //!
   TBranch        *b_ids;   //!
   TBranch        *b_p4_;   //!
   TBranch        *b_p4_fCoordinates_fPt;   //!
   TBranch        *b_p4_fCoordinates_fEta;   //!
   TBranch        *b_p4_fCoordinates_fPhi;   //!
   TBranch        *b_p4_fCoordinates_fE;   //!
   TBranch        *b_chargedHadronIsoR04_puppiNoMuonWeighted;   //!
   TBranch        *b_relativeIsoR03_withEA;   //!
   TBranch        *b_photonIsoR04_pfWeighted;   //!
   TBranch        *b_gen_p4_;   //!
   TBranch        *b_gen_p4_fCoordinates_fPt;   //!
   TBranch        *b_gen_p4_fCoordinates_fEta;   //!
   TBranch        *b_gen_p4_fCoordinates_fPhi;   //!
   TBranch        *b_gen_p4_fCoordinates_fE;   //!
   TBranch        *b_has_gen_particle;   //!
   TBranch        *b_chargedHadronIsoR03;   //!
      TBranch        *b_charge;   //! // test for signed char
   TBranch        *b_photonIsoR04_puppiWeighted;   //!
   TBranch        *b_y;   //!
   TBranch        *b_isLoose;   //!
   TBranch        *b_neutralHadronIsoR04_pfWeighted;   //!
   TBranch        *b_gen_y;   //!
   TBranch        *b_photonIsoR04;   //!
   TBranch        *b_puChargedHadronIsoR03;   //!
   TBranch        *b_isMedium;   //!
   TBranch        *b_relativeIsoR03_deltaBeta;   //!
   TBranch        *b_relativeIsoR03;   //!
   TBranch        *b_chargedHadronIsoR04;   //!
   TBranch        *b_neutralHadronIsoR04;   //!
   TBranch        *b_relativeIsoR04_deltaBeta;   //!
   TBranch        *b_puChargedHadronIsoR04;   //!
   TBranch        *b_relativeIsoR04;   //!
   TBranch        *b_relativeIsoR04_withEA;   //!

   base_muons(TTree *tree=0);
   virtual ~base_muons();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef base_muons_cxx
base_muons::base_muons(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150709/150709_160648/0000/output_mc_9.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150709/150709_160648/0000/output_mc_9.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150709/150709_160648/0000/output_mc_9.root:/muons");
      dir->GetObject("t",tree);

   }
   Init(tree);
}

base_muons::~base_muons()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t base_muons::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t base_muons::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void base_muons::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   relativeIsoR04_puppiNoMuonWeighted = 0;
   photonIsoR04_puppiNoMuonWeighted = 0;
   neutralHadronIsoR04_puppiNoMuonWeighted = 0;
   relativeIsoR04_puppiWeighted = 0;
   chargedHadronIsoR04_puppiWeighted = 0;
   relativeIsoR04_pfWeighted = 0;
   isTight = 0;
   neutralHadronIsoR04_puppiWeighted = 0;
   isHighPt = 0;
   isSoft = 0;
   photonIsoR03 = 0;
   neutralHadronIsoR03 = 0;
   //   gen_charge = 0;
   ids = 0;
   chargedHadronIsoR04_puppiNoMuonWeighted = 0;
   relativeIsoR03_withEA = 0;
   photonIsoR04_pfWeighted = 0;
   has_gen_particle = 0;
   chargedHadronIsoR03 = 0;
      charge = 0; // test for signed char
   photonIsoR04_puppiWeighted = 0;
   y = 0;
   isLoose = 0;
   neutralHadronIsoR04_pfWeighted = 0;
   gen_y = 0;
   photonIsoR04 = 0;
   puChargedHadronIsoR03 = 0;
   isMedium = 0;
   relativeIsoR03_deltaBeta = 0;
   relativeIsoR03 = 0;
   chargedHadronIsoR04 = 0;
   neutralHadronIsoR04 = 0;
   relativeIsoR04_deltaBeta = 0;
   puChargedHadronIsoR04 = 0;
   relativeIsoR04 = 0;
   relativeIsoR04_withEA = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("relativeIsoR04_puppiNoMuonWeighted", &relativeIsoR04_puppiNoMuonWeighted, &b_relativeIsoR04_puppiNoMuonWeighted);
   fChain->SetBranchAddress("photonIsoR04_puppiNoMuonWeighted", &photonIsoR04_puppiNoMuonWeighted, &b_photonIsoR04_puppiNoMuonWeighted);
   fChain->SetBranchAddress("neutralHadronIsoR04_puppiNoMuonWeighted", &neutralHadronIsoR04_puppiNoMuonWeighted, &b_neutralHadronIsoR04_puppiNoMuonWeighted);
   fChain->SetBranchAddress("relativeIsoR04_puppiWeighted", &relativeIsoR04_puppiWeighted, &b_relativeIsoR04_puppiWeighted);
   fChain->SetBranchAddress("chargedHadronIsoR04_puppiWeighted", &chargedHadronIsoR04_puppiWeighted, &b_chargedHadronIsoR04_puppiWeighted);
   fChain->SetBranchAddress("relativeIsoR04_pfWeighted", &relativeIsoR04_pfWeighted, &b_relativeIsoR04_pfWeighted);
   fChain->SetBranchAddress("isTight", &isTight, &b_isTight);
   fChain->SetBranchAddress("neutralHadronIsoR04_puppiWeighted", &neutralHadronIsoR04_puppiWeighted, &b_neutralHadronIsoR04_puppiWeighted);
   fChain->SetBranchAddress("isHighPt", &isHighPt, &b_isHighPt);
   fChain->SetBranchAddress("isSoft", &isSoft, &b_isSoft);
   fChain->SetBranchAddress("photonIsoR03", &photonIsoR03, &b_photonIsoR03);
   fChain->SetBranchAddress("neutralHadronIsoR03", &neutralHadronIsoR03, &b_neutralHadronIsoR03);
   //   fChain->SetBranchAddress("gen_charge", &gen_charge, &b_gen_charge);
   fChain->SetBranchAddress("ids", &ids, &b_ids);
   fChain->SetBranchAddress("p4", &p4_, &b_p4_);
   fChain->SetBranchAddress("p4.fCoordinates.fPt", p4_fCoordinates_fPt, &b_p4_fCoordinates_fPt);
   fChain->SetBranchAddress("p4.fCoordinates.fEta", p4_fCoordinates_fEta, &b_p4_fCoordinates_fEta);
   fChain->SetBranchAddress("p4.fCoordinates.fPhi", p4_fCoordinates_fPhi, &b_p4_fCoordinates_fPhi);
   fChain->SetBranchAddress("p4.fCoordinates.fE", p4_fCoordinates_fE, &b_p4_fCoordinates_fE);
   fChain->SetBranchAddress("chargedHadronIsoR04_puppiNoMuonWeighted", &chargedHadronIsoR04_puppiNoMuonWeighted, &b_chargedHadronIsoR04_puppiNoMuonWeighted);
   fChain->SetBranchAddress("relativeIsoR03_withEA", &relativeIsoR03_withEA, &b_relativeIsoR03_withEA);
   fChain->SetBranchAddress("photonIsoR04_pfWeighted", &photonIsoR04_pfWeighted, &b_photonIsoR04_pfWeighted);
   fChain->SetBranchAddress("gen_p4", &gen_p4_, &b_gen_p4_);
   fChain->SetBranchAddress("gen_p4.fCoordinates.fPt", gen_p4_fCoordinates_fPt, &b_gen_p4_fCoordinates_fPt);
   fChain->SetBranchAddress("gen_p4.fCoordinates.fEta", gen_p4_fCoordinates_fEta, &b_gen_p4_fCoordinates_fEta);
   fChain->SetBranchAddress("gen_p4.fCoordinates.fPhi", gen_p4_fCoordinates_fPhi, &b_gen_p4_fCoordinates_fPhi);
   fChain->SetBranchAddress("gen_p4.fCoordinates.fE", gen_p4_fCoordinates_fE, &b_gen_p4_fCoordinates_fE);
   fChain->SetBranchAddress("has_gen_particle", &has_gen_particle, &b_has_gen_particle);
   fChain->SetBranchAddress("chargedHadronIsoR03", &chargedHadronIsoR03, &b_chargedHadronIsoR03);
      fChain->SetBranchAddress("charge", &charge, &b_charge); // test for signed char
   fChain->SetBranchAddress("photonIsoR04_puppiWeighted", &photonIsoR04_puppiWeighted, &b_photonIsoR04_puppiWeighted);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("isLoose", &isLoose, &b_isLoose);
   fChain->SetBranchAddress("neutralHadronIsoR04_pfWeighted", &neutralHadronIsoR04_pfWeighted, &b_neutralHadronIsoR04_pfWeighted);
   fChain->SetBranchAddress("gen_y", &gen_y, &b_gen_y);
   fChain->SetBranchAddress("photonIsoR04", &photonIsoR04, &b_photonIsoR04);
   fChain->SetBranchAddress("puChargedHadronIsoR03", &puChargedHadronIsoR03, &b_puChargedHadronIsoR03);
   fChain->SetBranchAddress("isMedium", &isMedium, &b_isMedium);
   fChain->SetBranchAddress("relativeIsoR03_deltaBeta", &relativeIsoR03_deltaBeta, &b_relativeIsoR03_deltaBeta);
   fChain->SetBranchAddress("relativeIsoR03", &relativeIsoR03, &b_relativeIsoR03);
   fChain->SetBranchAddress("chargedHadronIsoR04", &chargedHadronIsoR04, &b_chargedHadronIsoR04);
   fChain->SetBranchAddress("neutralHadronIsoR04", &neutralHadronIsoR04, &b_neutralHadronIsoR04);
   fChain->SetBranchAddress("relativeIsoR04_deltaBeta", &relativeIsoR04_deltaBeta, &b_relativeIsoR04_deltaBeta);
   fChain->SetBranchAddress("puChargedHadronIsoR04", &puChargedHadronIsoR04, &b_puChargedHadronIsoR04);
   fChain->SetBranchAddress("relativeIsoR04", &relativeIsoR04, &b_relativeIsoR04);
   fChain->SetBranchAddress("relativeIsoR04_withEA", &relativeIsoR04_withEA, &b_relativeIsoR04_withEA);
   Notify();
}

Bool_t base_muons::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void base_muons::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t base_muons::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef base_muons_cxx

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 19 08:52:33 2016 by ROOT version 6.02/13
// from TTree t/t
// found on file: /eos/uscms/store/group/lpcmbja/jingli/PUID_76XMC/match_with_lepton/test.root
//////////////////////////////////////////////////////////

#ifndef PlotInput_h
#define PlotInput_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "Math/GenVector/PtEtaPhiE4D.h"
#include "vector"

const Int_t kMaxp4 = 1000;

class PlotInput {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<bool>    *allGenJet_patJetWithJetIDmatched;
   vector<bool>    *allGenJet_patJetmatched;
   vector<float>   *allGenJet_phi;
   vector<float>   *jetR;
   vector<float>   *jetRchg;
   vector<float>   *axisMinor;
   vector<float>   *axisMajor;
   vector<float>   *dRMatch;
   vector<int>     *nTot;
   vector<float>   *betaStar;
   vector<float>   *partonFlavor;
   vector<float>   *fRing3;
   vector<float>   *allGenJet_pt;
   vector<int>     *refpdgid;
   vector<float>   *fRing2;
   vector<float>   *DRweighted;
   vector<float>   *refarea;
   Int_t           p4_;
   Float_t         p4_fCoordinates_fPt[kMaxp4];   //[p4_]
   Float_t         p4_fCoordinates_fEta[kMaxp4];   //[p4_]
   Float_t         p4_fCoordinates_fPhi[kMaxp4];   //[p4_]
   Float_t         p4_fCoordinates_fE[kMaxp4];   //[p4_]
   vector<float>   *refdrjt;
   vector<float>   *pull;
   vector<unsigned short> *nCh;
   vector<float>   *fRing0;
   vector<float>   *allGenJet_m;
   vector<float>   *allGenJet_eta;
   vector<float>   *ptD;
   vector<float>   *beta;
   vector<float>   *fRing1;
   Float_t         scaleFactor;
   Float_t         rho;
   Int_t           nTrueInt;
   Int_t           npv;
   vector<float>   *jtpt;
   vector<float>   *jteta;
   vector<float>   *jtphi;
   vector<float>   *jtE;
   vector<float>   *refjtpt;
   vector<float>   *refjteta;
   vector<float>   *refjtphi;
   vector<float>   *refjtE;
   vector<float>   *dR_Jet_refJet;
   Float_t         weight;
   vector<float>   *dR_Jet_Muon;
   Int_t           n_TightMuon;
   Float_t         Muon1_charge;
   Float_t         Muon1_pt;
   Float_t         Muon1_eta;
   Float_t         Muon1_phi;
   Float_t         Muon1_E;
   Float_t         Muon2_charge;
   Float_t         Muon2_pt;
   Float_t         Muon2_eta;
   Float_t         Muon2_phi;
   Float_t         Muon2_E;
   Double_t        Mll;
   Int_t           b_atLeastOneMuon_TrigMatch;
   Float_t         DRTrigMatch1;
   Float_t         DRTrigMatch2;

   // List of branches
   TBranch        *b_allGenJet_patJetWithJetIDmatched;   //!
   TBranch        *b_allGenJet_patJetmatched;   //!
   TBranch        *b_allGenJet_phi;   //!
   TBranch        *b_jetR;   //!
   TBranch        *b_jetRchg;   //!
   TBranch        *b_axisMinor;   //!
   TBranch        *b_axisMajor;   //!
   TBranch        *b_dRMatch;   //!
   TBranch        *b_nTot;   //!
   TBranch        *b_betaStar;   //!
   TBranch        *b_partonFlavor;   //!
   TBranch        *b_fRing3;   //!
   TBranch        *b_allGenJet_pt;   //!
   TBranch        *b_refpdgid;   //!
   TBranch        *b_fRing2;   //!
   TBranch        *b_DRweighted;   //!
   TBranch        *b_refarea;   //!
   TBranch        *b_p4_;   //!
   TBranch        *b_p4_fCoordinates_fPt;   //!
   TBranch        *b_p4_fCoordinates_fEta;   //!
   TBranch        *b_p4_fCoordinates_fPhi;   //!
   TBranch        *b_p4_fCoordinates_fE;   //!
   TBranch        *b_refdrjt;   //!
   TBranch        *b_pull;   //!
   TBranch        *b_nCh;   //!
   TBranch        *b_fRing0;   //!
   TBranch        *b_allGenJet_m;   //!
   TBranch        *b_allGenJet_eta;   //!
   TBranch        *b_ptD;   //!
   TBranch        *b_beta;   //!
   TBranch        *b_fRing1;   //!
   TBranch        *b_scaleFactor;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtE;   //!
   TBranch        *b_refjtpt;   //!
   TBranch        *b_refjteta;   //!
   TBranch        *b_refjtphi;   //!
   TBranch        *b_refjtE;   //!
   TBranch        *b_dR_Jet_refJet;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_dR_Jet_Muon;   //!
   TBranch        *b_n_TightMuon;   //!
   TBranch        *b_Muon1_charge;   //!
   TBranch        *b_Muon1_pt;   //!
   TBranch        *b_Muon1_eta;   //!
   TBranch        *b_Muon1_phi;   //!
   TBranch        *b_Muon1_E;   //!
   TBranch        *b_Muon2_charge;   //!
   TBranch        *b_Muon2_pt;   //!
   TBranch        *b_Muon2_eta;   //!
   TBranch        *b_Muon2_phi;   //!
   TBranch        *b_Muon2_E;   //!
   TBranch        *b_Mll;   //!
   TBranch        *b_b_atLeastOneMuon_TrigMatch;   //!
   TBranch        *b_DRTrigMatch1;   //!
   TBranch        *b_DRTrigMatch2;   //!

   PlotInput(TTree *tree=0);
   virtual ~PlotInput();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PlotInput_cxx
PlotInput::PlotInput(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/uscms/store/group/lpcmbja/jingli/PUID_76XMC/match_with_lepton/test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/uscms/store/group/lpcmbja/jingli/PUID_76XMC/match_with_lepton/test.root");
      }
      f->GetObject("AK4PFCHS",tree);

   }
   Init(tree);
}

PlotInput::~PlotInput()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PlotInput::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PlotInput::LoadTree(Long64_t entry)
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

void PlotInput::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   allGenJet_patJetWithJetIDmatched = 0;
   allGenJet_patJetmatched = 0;
   allGenJet_phi = 0;
   jetR = 0;
   jetRchg = 0;
   axisMinor = 0;
   axisMajor = 0;
   dRMatch = 0;
   nTot = 0;
   betaStar = 0;
   partonFlavor = 0;
   fRing3 = 0;
   allGenJet_pt = 0;
   refpdgid = 0;
   fRing2 = 0;
   DRweighted = 0;
   refarea = 0;
   refdrjt = 0;
   pull = 0;
   nCh = 0;
   fRing0 = 0;
   allGenJet_m = 0;
   allGenJet_eta = 0;
   ptD = 0;
   beta = 0;
   fRing1 = 0;
   jtpt = 0;
   jteta = 0;
   jtphi = 0;
   jtE = 0;
   refjtpt = 0;
   refjteta = 0;
   refjtphi = 0;
   refjtE = 0;
   dR_Jet_refJet = 0;
   dR_Jet_Muon = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("allGenJet_patJetWithJetIDmatched", &allGenJet_patJetWithJetIDmatched, &b_allGenJet_patJetWithJetIDmatched);
   fChain->SetBranchAddress("allGenJet_patJetmatched", &allGenJet_patJetmatched, &b_allGenJet_patJetmatched);
   fChain->SetBranchAddress("allGenJet_phi", &allGenJet_phi, &b_allGenJet_phi);
   fChain->SetBranchAddress("jetR", &jetR, &b_jetR);
   fChain->SetBranchAddress("jetRchg", &jetRchg, &b_jetRchg);
   fChain->SetBranchAddress("axisMinor", &axisMinor, &b_axisMinor);
   fChain->SetBranchAddress("axisMajor", &axisMajor, &b_axisMajor);
   fChain->SetBranchAddress("dRMatch", &dRMatch, &b_dRMatch);
   fChain->SetBranchAddress("nTot", &nTot, &b_nTot);
   fChain->SetBranchAddress("betaStar", &betaStar, &b_betaStar);
   fChain->SetBranchAddress("partonFlavor", &partonFlavor, &b_partonFlavor);
   fChain->SetBranchAddress("fRing3", &fRing3, &b_fRing3);
   fChain->SetBranchAddress("allGenJet_pt", &allGenJet_pt, &b_allGenJet_pt);
   fChain->SetBranchAddress("refpdgid", &refpdgid, &b_refpdgid);
   fChain->SetBranchAddress("fRing2", &fRing2, &b_fRing2);
   fChain->SetBranchAddress("DRweighted", &DRweighted, &b_DRweighted);
   fChain->SetBranchAddress("refarea", &refarea, &b_refarea);
   fChain->SetBranchAddress("p4", &p4_, &b_p4_);
   fChain->SetBranchAddress("p4.fCoordinates.fPt", p4_fCoordinates_fPt, &b_p4_fCoordinates_fPt);
   fChain->SetBranchAddress("p4.fCoordinates.fEta", p4_fCoordinates_fEta, &b_p4_fCoordinates_fEta);
   fChain->SetBranchAddress("p4.fCoordinates.fPhi", p4_fCoordinates_fPhi, &b_p4_fCoordinates_fPhi);
   fChain->SetBranchAddress("p4.fCoordinates.fE", p4_fCoordinates_fE, &b_p4_fCoordinates_fE);
   fChain->SetBranchAddress("refdrjt", &refdrjt, &b_refdrjt);
   fChain->SetBranchAddress("pull", &pull, &b_pull);
   fChain->SetBranchAddress("nCh", &nCh, &b_nCh);
   fChain->SetBranchAddress("fRing0", &fRing0, &b_fRing0);
   fChain->SetBranchAddress("allGenJet_m", &allGenJet_m, &b_allGenJet_m);
   fChain->SetBranchAddress("allGenJet_eta", &allGenJet_eta, &b_allGenJet_eta);
   fChain->SetBranchAddress("ptD", &ptD, &b_ptD);
   fChain->SetBranchAddress("beta", &beta, &b_beta);
   fChain->SetBranchAddress("fRing1", &fRing1, &b_fRing1);
   fChain->SetBranchAddress("scaleFactor", &scaleFactor, &b_scaleFactor);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("jtpt", &jtpt, &b_jtpt);
   fChain->SetBranchAddress("jteta", &jteta, &b_jteta);
   fChain->SetBranchAddress("jtphi", &jtphi, &b_jtphi);
   fChain->SetBranchAddress("jtE", &jtE, &b_jtE);
   fChain->SetBranchAddress("refjtpt", &refjtpt, &b_refjtpt);
   fChain->SetBranchAddress("refjteta", &refjteta, &b_refjteta);
   fChain->SetBranchAddress("refjtphi", &refjtphi, &b_refjtphi);
   fChain->SetBranchAddress("refjtE", &refjtE, &b_refjtE);
   fChain->SetBranchAddress("dR_Jet_refJet", &dR_Jet_refJet, &b_dR_Jet_refJet);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("dR_Jet_Muon", &dR_Jet_Muon, &b_dR_Jet_Muon);
   fChain->SetBranchAddress("n_TightMuon", &n_TightMuon, &b_n_TightMuon);
   fChain->SetBranchAddress("Muon1_charge", &Muon1_charge, &b_Muon1_charge);
   fChain->SetBranchAddress("Muon1_pt", &Muon1_pt, &b_Muon1_pt);
   fChain->SetBranchAddress("Muon1_eta", &Muon1_eta, &b_Muon1_eta);
   fChain->SetBranchAddress("Muon1_phi", &Muon1_phi, &b_Muon1_phi);
   fChain->SetBranchAddress("Muon1_E", &Muon1_E, &b_Muon1_E);
   fChain->SetBranchAddress("Muon2_charge", &Muon2_charge, &b_Muon2_charge);
   fChain->SetBranchAddress("Muon2_pt", &Muon2_pt, &b_Muon2_pt);
   fChain->SetBranchAddress("Muon2_eta", &Muon2_eta, &b_Muon2_eta);
   fChain->SetBranchAddress("Muon2_phi", &Muon2_phi, &b_Muon2_phi);
   fChain->SetBranchAddress("Muon2_E", &Muon2_E, &b_Muon2_E);
   fChain->SetBranchAddress("Mll", &Mll, &b_Mll);
   fChain->SetBranchAddress("b_atLeastOneMuon_TrigMatch", &b_atLeastOneMuon_TrigMatch, &b_b_atLeastOneMuon_TrigMatch);
   fChain->SetBranchAddress("DRTrigMatch1", &DRTrigMatch1, &b_DRTrigMatch1);
   fChain->SetBranchAddress("DRTrigMatch2", &DRTrigMatch2, &b_DRTrigMatch2);
   Notify();
}

Bool_t PlotInput::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PlotInput::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PlotInput::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PlotInput_cxx

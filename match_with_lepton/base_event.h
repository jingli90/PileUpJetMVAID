//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 17 11:12:50 2015 by ROOT version 6.02/05
// from TTree t/t
// found on file: /eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150709/150709_160648/0000/output_mc_9.root
//////////////////////////////////////////////////////////

#ifndef base_event_h
#define base_event_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

#include "base.h" // put the size in to one file to avoid duplication
				  // It seems that they are not used in event/t

class base_event {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nMEPartonsFiltered;
 //pair<float,float> *pdf_x;
   Float_t         first_pdf;
   Float_t         second_pdf;
 //pair<int,int>   *pdf_id;
   Int_t           first_id;
   Int_t           second_id;
   std::vector<int>     *bunch_spacings;
   ULong64_t       evt;
   Float_t         alphaQCD;
   Int_t           nMEPartons;
   ULong64_t       run;
   Float_t         alphaQED;
   std::vector<float>   *rhos;
   ULong64_t       lumi;
   std::vector<float>   *pu_sumpt_lowpt;
   std::vector<int>     *npus;
   ULong64_t       npv;
   Float_t         rho;
   Float_t         pthat;
   Float_t         qScale;
   Int_t           npuOOT;
   std::vector<float>   *tnpus;
   Float_t         weight;
   std::vector<int>     *pu_ntrks_lowpt;
   std::vector<float>   *pu_sumpt_highpt;
   std::vector<std::vector<float> > *pu_zpositions;
   std::vector<int>     *pu_ntrks_highpt;
   Int_t           npuIT;
   std::vector<int>     *bxns;
   Int_t           nTrueInt;

   // List of branches
   TBranch        *b_nMEPartonsFiltered;   //!
   TBranch        *b_pdf_x_first;   //!
   TBranch        *b_pdf_x_second;   //!
   TBranch        *b_pdf_id_first;   //!
   TBranch        *b_pdf_id_second;   //!
   TBranch        *b_bunch_spacings;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_alphaQCD;   //!
   TBranch        *b_nMEPartons;   //!
   TBranch        *b_run;   //!
   TBranch        *b_alphaQED;   //!
   TBranch        *b_rhos;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_pu_sumpt_lowpt;   //!
   TBranch        *b_npus;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_qScale;   //!
   TBranch        *b_npuOOT;   //!
   TBranch        *b_tnpus;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_pu_ntrks_lowpt;   //!
   TBranch        *b_pu_sumpt_highpt;   //!
   TBranch        *b_pu_zpositions;   //!
   TBranch        *b_pu_ntrks_highpt;   //!
   TBranch        *b_npuIT;   //!
   TBranch        *b_bxns;   //!
   TBranch        *b_nTrueInt;   //!

   base_event(TTree *tree=0);
   virtual ~base_event();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef base_event_cxx
base_event::base_event(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150709/150709_160648/0000/output_mc_9.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150709/150709_160648/0000/output_mc_9.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150709/150709_160648/0000/output_mc_9.root:/event");
      dir->GetObject("t",tree);

   }
   Init(tree);
}

base_event::~base_event()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t base_event::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t base_event::LoadTree(Long64_t entry)
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

void base_event::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   bunch_spacings = 0;
   rhos = 0;
   pu_sumpt_lowpt = 0;
   npus = 0;
   tnpus = 0;
   pu_ntrks_lowpt = 0;
   pu_sumpt_highpt = 0;
   pu_zpositions = 0;
   pu_ntrks_highpt = 0;
   bxns = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nMEPartonsFiltered", &nMEPartonsFiltered, &b_nMEPartonsFiltered);
   fChain->SetBranchAddress("first", &first_pdf, &b_pdf_x_first);
   fChain->SetBranchAddress("second", &second_pdf, &b_pdf_x_second);
   fChain->SetBranchAddress("first", &first_id, &b_pdf_id_first);
   fChain->SetBranchAddress("second", &second_id, &b_pdf_id_second);
   fChain->SetBranchAddress("bunch_spacings", &bunch_spacings, &b_bunch_spacings);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("alphaQCD", &alphaQCD, &b_alphaQCD);
   fChain->SetBranchAddress("nMEPartons", &nMEPartons, &b_nMEPartons);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("alphaQED", &alphaQED, &b_alphaQED);
   fChain->SetBranchAddress("rhos", &rhos, &b_rhos);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("pu_sumpt_lowpt", &pu_sumpt_lowpt, &b_pu_sumpt_lowpt);
   fChain->SetBranchAddress("npus", &npus, &b_npus);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
   fChain->SetBranchAddress("npuOOT", &npuOOT, &b_npuOOT);
   fChain->SetBranchAddress("tnpus", &tnpus, &b_tnpus);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("pu_ntrks_lowpt", &pu_ntrks_lowpt, &b_pu_ntrks_lowpt);
   fChain->SetBranchAddress("pu_sumpt_highpt", &pu_sumpt_highpt, &b_pu_sumpt_highpt);
   fChain->SetBranchAddress("pu_zpositions", &pu_zpositions, &b_pu_zpositions);
   fChain->SetBranchAddress("pu_ntrks_highpt", &pu_ntrks_highpt, &b_pu_ntrks_highpt);
   fChain->SetBranchAddress("npuIT", &npuIT, &b_npuIT);
   fChain->SetBranchAddress("bxns", &bxns, &b_bxns);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   Notify();
}

Bool_t base_event::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void base_event::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t base_event::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef base_event_cxx

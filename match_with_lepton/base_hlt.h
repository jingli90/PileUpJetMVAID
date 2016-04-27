//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 27 15:34:44 2015 by ROOT version 6.02/05
// from TTree t/t
// found on file: /eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150722-v2_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15DR74-StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/150722_223310/0000/output_mc_1.root
//////////////////////////////////////////////////////////

#ifndef base_hlt_h
#define base_hlt_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class base_hlt {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   std::vector<std::vector<bool> > *objects_paths_islast;
   std::vector<std::vector<int> > *objects_paths;
   std::vector<std::vector<bool> > *objects_paths_isl3;
   std::vector<float>   *objects_phi;
   std::vector<float>   *objects_e;
   std::vector<float>   *objects_eta;
   std::vector<float>   *objects_pt;
   std::vector<int>     *prescales;
   std::vector<std::string>  *paths;

   // List of branches
   TBranch        *b_objects_paths_islast;   //!
   TBranch        *b_objects_paths;   //!
   TBranch        *b_objects_paths_isl3;   //!
   TBranch        *b_objects_phi;   //!
   TBranch        *b_objects_e;   //!
   TBranch        *b_objects_eta;   //!
   TBranch        *b_objects_pt;   //!
   TBranch        *b_prescales;   //!
   TBranch        *b_paths;   //!

   base_hlt(TTree *tree=0);
   virtual ~base_hlt();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef base_hlt_cxx
base_hlt::base_hlt(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150722-v2_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15DR74-StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/150722_223310/0000/output_mc_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150722-v2_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15DR74-StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/150722_223310/0000/output_mc_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/eos/uscms/store/group/lpcmbja/jingli/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator746p2_150722-v2_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15DR74-StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/150722_223310/0000/output_mc_1.root:/hlt");
      dir->GetObject("t",tree);

   }
   Init(tree);
}

base_hlt::~base_hlt()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t base_hlt::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t base_hlt::LoadTree(Long64_t entry)
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

void base_hlt::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   objects_paths_islast = 0;
   objects_paths = 0;
   objects_paths_isl3 = 0;
   objects_phi = 0;
   objects_e = 0;
   objects_eta = 0;
   objects_pt = 0;
   prescales = 0;
   paths = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("objects_paths_islast", &objects_paths_islast, &b_objects_paths_islast);
   fChain->SetBranchAddress("objects_paths", &objects_paths, &b_objects_paths);
   fChain->SetBranchAddress("objects_paths_isl3", &objects_paths_isl3, &b_objects_paths_isl3);
   fChain->SetBranchAddress("objects_phi", &objects_phi, &b_objects_phi);
   fChain->SetBranchAddress("objects_e", &objects_e, &b_objects_e);
   fChain->SetBranchAddress("objects_eta", &objects_eta, &b_objects_eta);
   fChain->SetBranchAddress("objects_pt", &objects_pt, &b_objects_pt);
   fChain->SetBranchAddress("prescales", &prescales, &b_prescales);
   fChain->SetBranchAddress("paths", &paths, &b_paths);
   Notify();
}

Bool_t base_hlt::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void base_hlt::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t base_hlt::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef base_hlt_cxx

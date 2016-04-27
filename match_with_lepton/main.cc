#include <iostream>
#include <string>

#include <vector>

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TROOT.h>

#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <fstream>
#include <math.h>

#include <TLorentzVector.h>

#include <iostream>
using namespace std;

#include "base.h"
#include "base_jet.h"
#include "base_event.h"
#include "base_hlt.h"
#include "base_muons.h"
#include "base_electrons.h"

double calc_dr2( double e1, double e2, double p1, double p2 ){
	//cout<<"M_PI = "<<M_PI<<endl;
	double d_phi = fabs( p1 - p2 );
	d_phi = d_phi < M_PI ? d_phi : 2*M_PI - d_phi ;
	double d_eta = e1 - e2 ;
	return d_eta * d_eta + d_phi * d_phi ;
}

int main(int argc,char *argv[]){  

  TString JetAlgo = "AK4PFCHS";
//  TString JetAlgo = "AK4PF";
  cout<<"JetAlgo="<<JetAlgo<<endl;

  Float_t scaleFactor;

  TString filename = "/eos/uscms/store/group/lpcmbja/jingli/PUID_76XMC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_jmevalidator763p2_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160419_093132/output_mc.root";
  TString outfilename = "/eos/uscms/store/group/lpcmbja/jingli/PUID_76XMC/match_with_lepton/test.root";
  scaleFactor=1.;
  
  TChain *c_jet  = new TChain(Form("%s/t",JetAlgo.Data()));
  TChain *c_event     = new TChain("event/t");
  TChain *c_hlt       = new TChain("hlt/t");
  TChain *c_muons     = new TChain("muons/t");
//  TChain *c_electrons = new TChain("electrons/t");

  c_jet->Add(filename.Data());
  c_event->Add(filename.Data());
  c_hlt->Add(filename.Data());
  c_muons->Add(filename.Data());
//  c_electrons->Add(filename.Data());

  if (c_jet == 0){
	  cout<<"error! c_jet not init";
	  return 1;
  }

  base_jet b_jet(c_jet);
  base_event b_event(c_event);
  base_hlt b_hlt(c_hlt);
  base_muons b_muons(c_muons);
//  base_electrons b_electrons(c_electrons);

  c_jet->SetBranchStatus("*",0);
  c_jet->SetBranchStatus("DRweighted",1);
  c_jet->SetBranchStatus("nTot",1);
  c_jet->SetBranchStatus("nCh",1);
  c_jet->SetBranchStatus("axisMajor",1);
  c_jet->SetBranchStatus("axisMinor",1);
  c_jet->SetBranchStatus("fRing0",1);
  c_jet->SetBranchStatus("fRing1",1);
  c_jet->SetBranchStatus("fRing2",1);
  c_jet->SetBranchStatus("fRing3",1);
  c_jet->SetBranchStatus("ptD",1);
  c_jet->SetBranchStatus("beta",1);
  c_jet->SetBranchStatus("betaStar",1);
  c_jet->SetBranchStatus("pull",1);
  c_jet->SetBranchStatus("jetR",1);
  c_jet->SetBranchStatus("jetRchg",1);
  c_jet->SetBranchStatus("p4*",1);
  c_jet->SetBranchStatus("dRMatch",1);
  c_jet->SetBranchStatus("allGenJet*",1);
  c_jet->SetBranchStatus("ref*",1);
  c_jet->SetBranchStatus("partonFlavor",1);

  //Long64_t nentries=0, nentries_jet=0, nentries_event=0, nentries_hlt=0, nentries_muons=0, nentries_electrons=0;
  Long64_t nentries=0, nentries_jet=0, nentries_event=0, nentries_hlt=0, nentries_muons=0;
  nentries_jet = c_jet->GetEntries();
  nentries_event = c_event->GetEntries();
  nentries_hlt=c_hlt->GetEntries();
  nentries_muons = c_muons->GetEntries();
//  nentries_electrons = c_electrons->GetEntries();
  
  nentries=nentries_jet;

  //if((nentries_jet!=nentries_event)||(nentries_jet!=nentries_hlt)||(nentries_jet!=nentries_muons)||(nentries_jet!=nentries_electrons))
  if((nentries_jet!=nentries_event)||(nentries_jet!=nentries_hlt)||(nentries_jet!=nentries_muons))
  {
	  cout<<"error! Three trees has different entries:"<<endl;
	  //cout<<"jet: "<<nentries_jet<<" event: "<<nentries_event<<" hlt: "<<nentries_hlt<<" muons: "<<nentries_muons<<" electrons: "<<nentries_electrons<<endl;
	  cout<<"jet: "<<nentries_jet<<" event: "<<nentries_event<<" hlt: "<<nentries_hlt<<" muons: "<<nentries_muons<<endl;
//	  return 1;
  }
  else{
	  cout<<"total: "<<nentries<<endl;
  }

  TFile * fout = new TFile(outfilename.Data(),"RECREATE");
  TTree * tout = c_jet->CloneTree(0);
 
  tout->Branch("scaleFactor",&scaleFactor,"scaleFactor/F");
  Float_t         rho;
  tout->Branch("rho",&rho,"rho/F");
  Int_t           nTrueInt;
  tout->Branch("nTrueInt",&nTrueInt,"nTrueInt/I");
  Int_t npv;
  tout->Branch("npv",&npv,"npv/I");

  std::vector<float>* jtpt;
  std::vector<float>* jteta;
  std::vector<float>* jtphi;
  std::vector<float>* jtE;
  std::vector<float>* refjtpt;
  std::vector<float>* refjteta;
  std::vector<float>* refjtphi;
  std::vector<float>* refjtE;
  jtpt = new std::vector<float>;
  jteta = new std::vector<float>;
  jtphi = new std::vector<float>;
  jtE = new std::vector<float>;
  refjtpt = new std::vector<float>;
  refjteta = new std::vector<float>;
  refjtphi = new std::vector<float>;
  refjtE = new std::vector<float>;
  tout->Branch("jtpt","vector<Float_t>",&jtpt);
  tout->Branch("jteta","vector<Float_t>",&jteta);
  tout->Branch("jtphi","vector<Float_t>",&jtphi);
  tout->Branch("jtE","vector<Float_t>",&jtE);
  tout->Branch("refjtpt","vector<Float_t>",&refjtpt);
  tout->Branch("refjteta","vector<Float_t>",&refjteta); 
  tout->Branch("refjtphi","vector<Float_t>",&refjtphi);
  tout->Branch("refjtE","vector<Float_t>",&refjtE);

  Float_t jtpt_,jteta_,jtphi_,jtE_;
  TLorentzVector vJet;

  Float_t refjtpt_,refjteta_,refjtphi_,refjtE_;
  TLorentzVector vrefJet;

  std::vector<float>* dR_Jet_refJet;
  dR_Jet_refJet = new std::vector<float>;
  tout->Branch("dR_Jet_refJet","vector<Float_t>",&dR_Jet_refJet);
  Float_t dR_Jet_refJet_;

  Float_t weight;
  tout->Branch("weight",&weight,"weight/F");

  std::vector<float>* dR_Jet_Muon;
  dR_Jet_Muon = new std::vector<float>;
  tout->Branch("dR_Jet_Muon","vector<Float_t>",&dR_Jet_Muon);
  Float_t dR_Jet_Muon_; 

  Float_t Muoncharge_;
  Int_t n_TightMuon=0;
  Float_t Muon1_charge;
  Float_t Muon1_pt, Muon1_eta, Muon1_phi, Muon1_E;
  Float_t Muon2_charge;
  Float_t Muon2_pt, Muon2_eta, Muon2_phi, Muon2_E;
  Double_t Mll;
  Int_t b_atLeastOneMuon_TrigMatch = 0 ; 
  Float_t DRTrigMatch1=9999,DRTrigMatch2=9999;
  Float_t DRTrigMatch[2]={9999,9999};

  tout->Branch("n_TightMuon",&n_TightMuon,"n_TightMuon/I");
  tout->Branch("Muon1_charge",&Muon1_charge,"Muon1_charge/F");
  tout->Branch("Muon1_pt",&Muon1_pt,"Muon1_pt/F");
  tout->Branch("Muon1_eta",&Muon1_eta,"Muon1_eta/F");
  tout->Branch("Muon1_phi",&Muon1_phi,"Muon1_phi/F");
  tout->Branch("Muon1_E",&Muon1_E,"Muon1_E/F");
  tout->Branch("Muon2_charge",&Muon2_charge,"Muon2_charge/F");
  tout->Branch("Muon2_pt",&Muon2_pt,"Muon2_pt/F");
  tout->Branch("Muon2_eta",&Muon2_eta,"Muon2_eta/F");
  tout->Branch("Muon2_phi",&Muon2_phi,"Muon2_phi/F");
  tout->Branch("Muon2_E",&Muon2_E,"Muon2_E/F");
  tout->Branch("Mll",&Mll,"Mll/D");
  tout->Branch("b_atLeastOneMuon_TrigMatch",&b_atLeastOneMuon_TrigMatch,"b_atLeastOneMuon_TrigMatch/I");
  tout->Branch("DRTrigMatch1",&DRTrigMatch1,"DRTrigMatch1/F");
  tout->Branch("DRTrigMatch2",&DRTrigMatch2,"DRTrigMatch2/F");

  TH1F h_Mll("Mll","Mll", 200 , 50 , 150 ) ;
  TH1F h_pt_jet_match_with_muon("pt_jet_match_with_muon","pt_jet_match_with_muon",60,0,120);
  TH1F h_pt_jet_NOTmatch_with_muon("pt_jet_NOTmatch_with_muon","pt_jet_NOTmatch_with_muon",60,0,120);
  TH1F h_pt_jet_Zmumu("pt_jet_Zmumu","pt_jet_Zmumu",60,0,120);

  TLorentzVector v0(0,0,0,0);

  const std::string requiredTriggerMenu("HLT_IsoMu20_v");

  for (Long64_t jentry=0; jentry<nentries;jentry++)
  //for (Long64_t jentry=0; jentry<10000;jentry++)
  {
	  if (jentry%1000==0)cout<<"entry: "<<jentry<<endl;
	  //cout<<"entry: "<<jentry<<endl;

	  jtpt->clear(); 
	  jteta->clear();  
	  jtphi->clear();  
	  jtE->clear(); 
	  refjtpt->clear(); 
	  refjteta->clear();
	  refjtphi->clear();
	  refjtE->clear();
	
	  dR_Jet_Muon->clear();
	  dR_Jet_refJet->clear();
	
	  for(unsigned int i=0;i<kMaxp4;i++){
	  	if(b_jet.p4_fCoordinates_fPt[i]!=0){
	  		b_jet.p4_fCoordinates_fPt[i]=0;
	  		b_jet.p4_fCoordinates_fEta[i]=0;
	  		b_jet.p4_fCoordinates_fPhi[i]=0;
	  		b_jet.p4_fCoordinates_fE[i]=0;
	  	}
	  }
	  for(unsigned int i=0;i<kMaxgen_p4;i++){
	  	if(b_jet.gen_p4_fCoordinates_fPt[i]!=0){
	  		b_jet.gen_p4_fCoordinates_fPt[i]=0;
	  		b_jet.gen_p4_fCoordinates_fEta[i]=0;
	  		b_jet.gen_p4_fCoordinates_fPhi[i]=0;
	  		b_jet.gen_p4_fCoordinates_fE[i]=0;
	  	}
	  }	  
	  for(unsigned int i=0;i<kMaxp4;i++){
	  	if(b_muons.p4_fCoordinates_fPt[i]!=0){
	  		b_muons.p4_fCoordinates_fPt[i]=0;
	  		b_muons.p4_fCoordinates_fEta[i]=0;
	  		b_muons.p4_fCoordinates_fPhi[i]=0;
	  		b_muons.p4_fCoordinates_fE[i]=0;
	  	}
	  }
	  for(unsigned int i=0;i<kMaxgen_p4;i++){
	  	if(b_muons.gen_p4_fCoordinates_fPt[i]!=0){
	  		b_muons.gen_p4_fCoordinates_fPt[i]=0;
	  		b_muons.gen_p4_fCoordinates_fEta[i]=0;
	  		b_muons.gen_p4_fCoordinates_fPhi[i]=0;
	  		b_muons.gen_p4_fCoordinates_fE[i]=0;
	  	}
	  }	  

	  c_jet->GetEntry(jentry);
	  c_event->GetEntry(jentry);
	  c_muons->GetEntry(jentry);
	  c_hlt->GetEntry(jentry);

	  //cout<<"first id: "<<b_event.first_id<<" second id: "<<b_event.second_id<<endl;
	  //cout<<"jet size: "<<b_jet.y->size()<<" "<<b_jet.p4_<<endl;
	  //cout<<"muon size: "<<b_muons.y->size()<<" "<<b_muons.p4_<<endl;
	  //cout<<"electrons size: "<<b_electrons.y->size()<<" "<<b_electrons.p4_<<endl;

	  std::vector< int > TrigObjIdx ; 
	  for( unsigned int i = 0 ; i < b_hlt.paths->size(); i ++){
		  std::string * trig = & b_hlt.paths->at(i);
		  if( trig->find( requiredTriggerMenu ) == std::string::npos ) continue ;
		  //cout<<"i = "<<i<<" hlt: "<<(*trig)<<endl;
		  for( int iObj = 0 ; iObj < b_hlt.objects_paths->size() ; iObj ++ ){
			  // Search obje that matches to the trigger.
			  //cout<<"iObj "<<iObj<<endl;
			  std::vector<int> p = b_hlt.objects_paths ->at(iObj) ;
			  std::vector<bool> b = b_hlt.objects_paths_islast->at(iObj) ;
			  for( unsigned int j = 0 ; j < p.size() ; j++ ){
				  //cout<<"j "<<j<<endl;
				  //cout<<"b_hlt.objects_paths->at("<<iObj<<").at("<<j<<") = "<<p.at(j)<<endl;
				  //cout<<"b_hlt.objects_paths_islast->at("<<iObj<<").at("<<j<<") = "<<b.at(j)<<endl;
				  if( b.at(j) && ( p.at(j) == i )) {
					  TrigObjIdx . push_back( iObj );
					  //cout<<"b_hlt.objects_paths->at("<<iObj<<").at("<<j<<") = "<<p.at(j)<<endl;
					  //cout<<"b_hlt.objects_paths_islast->at("<<iObj<<").at("<<j<<") = "<<b.at(j)<<endl;
				  }
			  }
		  }
	  }

	  if( TrigObjIdx .size() == 0 ) continue ;

	  // exact 2 tight muon
	  n_TightMuon = 0 ;
	  long idx_mu[2] = { -1 , -1 };
	  for (int i=0; i<b_muons.p4_; i++){
		  Muoncharge_ = b_muons.charge->at(i);
		  //cout<<Muoncharge_<<endl;

		  if(b_muons.isTight->at(i)){
			  if( idx_mu[ 0 ] == -1 ){
				  idx_mu[0] = i ;
			  }else{
				  idx_mu[1] = i ;
			  }
			  n_TightMuon ++ ; 
		  }
	  }
	  //cout<<n_TightMuon<<endl;
	  if( n_TightMuon != 2 ) continue ;

	  Muon1_charge=b_muons.charge->at(idx_mu[0]);
	  Muon2_charge=b_muons.charge->at(idx_mu[1]);

	  Muon1_pt=b_muons.p4_fCoordinates_fPt[idx_mu[0]];
	  Muon1_eta=b_muons.p4_fCoordinates_fEta[idx_mu[0]];
	  Muon1_phi=b_muons.p4_fCoordinates_fPhi[idx_mu[0]];
	  Muon1_E=b_muons.p4_fCoordinates_fE[idx_mu[0]];

	  Muon2_pt=b_muons.p4_fCoordinates_fPt[idx_mu[1]];
	  Muon2_eta=b_muons.p4_fCoordinates_fEta[idx_mu[1]];
	  Muon2_phi=b_muons.p4_fCoordinates_fPhi[idx_mu[1]];
	  Muon2_E=b_muons.p4_fCoordinates_fE[idx_mu[1]];

	  if(Muon1_charge*Muon2_charge>0){
		  std::cout <<"Same sign --- reject. " << std::endl ; 
		  continue;
	  }
	  TLorentzVector v_mu1 ; 
	  v_mu1.SetPtEtaPhiE(
				  b_muons.p4_fCoordinates_fPt[idx_mu[0]],
				  b_muons.p4_fCoordinates_fEta[idx_mu[0]],
				  b_muons.p4_fCoordinates_fPhi[idx_mu[0]],
				  b_muons.p4_fCoordinates_fE[idx_mu[0]] 
				  );
	  TLorentzVector v_mu2 ; 
	  v_mu2.SetPtEtaPhiE(
				  b_muons.p4_fCoordinates_fPt[idx_mu[1]],
				  b_muons.p4_fCoordinates_fEta[idx_mu[1]],
				  b_muons.p4_fCoordinates_fPhi[idx_mu[1]],
				  b_muons.p4_fCoordinates_fE[idx_mu[1]] 
				  );
	  Mll = ( v_mu1 + v_mu2 ).M() ; 
	  h_Mll . Fill( Mll) ;
	  //cout<<Mll<<endl;
	  //if( fabs( Mll - 91.0 ) > 10 ) continue ; // original
	  if(!(Mll>70 && Mll<110)) continue; //  Kostas 
	  b_atLeastOneMuon_TrigMatch=0;
	  const float ThresholdTrigMatchDR2 = 0.1 * 0.1 ; 
	  DRTrigMatch[0]=999;
	  DRTrigMatch[1]=999;
	  for( int i = 0 ; i < 2 ; i ++ ){
		  int idx_other = ( i== 0 ? 1 : 0 );
		  double trigmatch_dr2 = 999 * 999;
		  double trigmatch_dr=999;
		  bool trig_match = false ;
		  for( unsigned int iTrig = 0 ; iTrig < TrigObjIdx.size() ; iTrig ++ ){
			  trigmatch_dr2 = calc_dr2(
						  b_muons.p4_fCoordinates_fEta[idx_mu[idx_other]],
						  b_hlt. objects_eta ->at( TrigObjIdx.at( iTrig ) ),
						  b_muons.p4_fCoordinates_fPhi[idx_mu[idx_other]],
						  b_hlt. objects_phi ->at( TrigObjIdx.at( iTrig ) )
						  );
			  trigmatch_dr=sqrt(trigmatch_dr2);
			  if(trigmatch_dr<DRTrigMatch[idx_other]){
				  DRTrigMatch[idx_other]=trigmatch_dr;
			  }
			  if( trigmatch_dr2 < ThresholdTrigMatchDR2 ) trig_match = true ; 
		  }
		  if( ! trig_match ) continue ;
		  b_atLeastOneMuon_TrigMatch = true ; 
	  }
	  DRTrigMatch1=DRTrigMatch[0];
	  DRTrigMatch2=DRTrigMatch[1];

	  for (int i=0; i<b_jet.p4_; i++){
		  jtpt_= b_jet.p4_fCoordinates_fPt[i];
		  jteta_=b_jet.p4_fCoordinates_fEta[i];
		  jtphi_=b_jet.p4_fCoordinates_fPhi[i];
		  jtE_=  b_jet.p4_fCoordinates_fE[i];
		  vJet.SetPtEtaPhiE(jtpt_,jteta_,jtphi_,jtE_);

		  refjtpt_= b_jet.gen_p4_fCoordinates_fPt[i];
		  refjteta_=b_jet.gen_p4_fCoordinates_fEta[i];
		  refjtphi_=b_jet.gen_p4_fCoordinates_fPhi[i];
		  refjtE_=  b_jet.gen_p4_fCoordinates_fE[i];
		  vrefJet.SetPtEtaPhiE(refjtpt_,refjteta_,refjtphi_,refjtE_);
	
		  jtpt->push_back(jtpt_);
		  jteta->push_back(jteta_);
		  jtphi->push_back(jtphi_);
		  jtE->push_back(jtE_);
		  refjtpt->push_back(refjtpt_);
		  refjteta->push_back(refjteta_);
		  refjtphi->push_back(refjtphi_);
		  refjtE->push_back(refjtE_);

//		  dR_refJet_refMuon_=999;
		  dR_Jet_Muon_=999;

		  if(vJet!=v0 && vJet.DeltaR(v_mu1)<dR_Jet_Muon_){
			  dR_Jet_Muon_=vJet.DeltaR(v_mu1);
		  }

		  if(vJet!=v0 && vJet.DeltaR(v_mu2)<dR_Jet_Muon_){
			  dR_Jet_Muon_=vJet.DeltaR(v_mu2);
		  }

		  dR_Jet_Muon->push_back(dR_Jet_Muon_);
//		  //cout<<"dR_refJet_refMuon_="<<dR_refJet_refMuon_<<endl;
//		  cout<<"dR_Jet_Muon_="<<dR_Jet_Muon_<<endl;

		  if(b_atLeastOneMuon_TrigMatch && jtpt_>15){
			  h_pt_jet_Zmumu.Fill(jtpt_);
			  if(dR_Jet_Muon_<0.4) h_pt_jet_match_with_muon.Fill(jtpt_);
			  if(dR_Jet_Muon_>0.4) h_pt_jet_NOTmatch_with_muon.Fill(jtpt_);
		  }
	  }

	  for (int i=0; i<b_jet.p4_; i++){
		  dR_Jet_refJet_=999;
		  vJet.SetPtEtaPhiE(jtpt->at(i),jteta->at(i),jtphi->at(i),jtE->at(i));
		  for(int i_refjet=0;i_refjet<b_jet.p4_;i_refjet++){
			  vrefJet.SetPtEtaPhiE(refjtpt->at(i_refjet),refjteta->at(i_refjet),refjtphi->at(i_refjet),refjtE->at(i_refjet));
			  if(vJet!=v0 && vrefJet!=v0 && vJet.DeltaR(vrefJet)<dR_Jet_refJet_){
				  dR_Jet_refJet_=vJet.DeltaR(vrefJet);
				  //cout<<dR_Jet_refJet_<<endl;
			  }
		  }
		  dR_Jet_refJet->push_back(dR_Jet_refJet_);
	  }

	  rho=b_event.rho;
	  nTrueInt=b_event.nTrueInt;
	  npv=b_event.npv;

	  weight=b_event.weight;

	  tout->Fill();
  }

  tout->Write(JetAlgo.Data());
  h_Mll.Write();
  h_pt_jet_match_with_muon.Write();
  h_pt_jet_NOTmatch_with_muon.Write();
  h_pt_jet_Zmumu.Write();
  fout->Close();

  cout<<"Writing "<<JetAlgo<<" into "<<outfilename<<endl;

  return 0;
}

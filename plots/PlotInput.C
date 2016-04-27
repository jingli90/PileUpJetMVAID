#define PlotInput_cxx
#include "PlotInput.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
using namespace std;

void PlotInput::Loop()
{
   if (fChain == 0) return;

   TH1D * hDRweighted[3][3][4];
   TH1D * hnpv[3][3][4];
   TH1D * hnTot[3][3][4];
   TH1D * hnCh[3][3][4];
   TH1D * haxisMajor[3][3][4];
   TH1D * haxisMinor[3][3][4];
   TH1D * hfRing0[3][3][4];
   TH1D * hfRing1[3][3][4];
   TH1D * hfRing2[3][3][4];
   TH1D * hfRing3[3][3][4];
   TH1D * hptD[3][3][4];
   TH1D * hbeta[3][3][4];
   TH1D * hpull[3][3][4];
   TH1D * hjetR[3][3][4];
   TH1D * hjetRchg[3][3][4];
   TH1D * hjtpt[3][3][4];
   TH1D * hjteta[3][3][4];
   for(int i=0;i<3;i++){
	   for(int j=0;j<3;j++){
		   for(int k=0;k<4;k++){
			   hDRweighted[i][j][k]=new TH1D(Form("hDRweighted_%d_%d_%d",i,j,k),Form("hDRweighted_%d_%d_%d",i,j,k),100,0.,0.1);
			   hDRweighted[i][j][k]->SetTitle(";dR2Mean;");
			   hnpv[i][j][k]=new TH1D(Form("hnpv_%d_%d_%d",i,j,k),Form("hnpv_%d_%d_%d",i,j,k),100,0.,100.);
			   hnpv[i][j][k]->SetTitle(";npv;");
			   hnTot[i][j][k]=new TH1D(Form("hnTot_%d_%d_%d",i,j,k),Form("hnTot_%d_%d_%d",i,j,k),60,0.,60.);
			   hnTot[i][j][k]->SetTitle(";nParticles;");
			   hnCh[i][j][k]=new TH1D(Form("hnCh_%d_%d_%d",i,j,k),Form("hnCh_%d_%d_%d",i,j,k),40,0.,40.);
			   hnCh[i][j][k]->SetTitle(";nCharged;");
			   haxisMajor[i][j][k]=new TH1D(Form("haxisMajor_%d_%d_%d",i,j,k),Form("haxisMajor_%d_%d_%d",i,j,k),60,0.,0.3);
			   haxisMajor[i][j][k]->SetTitle(";majW;");
			   haxisMinor[i][j][k]=new TH1D(Form("haxisMinor_%d_%d_%d",i,j,k),Form("haxisMinor_%d_%d_%d",i,j,k),40,0.,0.2);
			   haxisMinor[i][j][k]->SetTitle(";minW;");
			   hfRing0[i][j][k]=new TH1D(Form("hfRing0_%d_%d_%d",i,j,k),Form("hfRing0_%d_%d_%d",i,j,k),100,0.,1.);
			   hfRing0[i][j][k]->SetTitle(";frac01;");
			   hfRing1[i][j][k]=new TH1D(Form("hfRing1_%d_%d_%d",i,j,k),Form("hfRing1_%d_%d_%d",i,j,k),100,0.,1.);
			   hfRing1[i][j][k]->SetTitle(";frac02;");
			   hfRing2[i][j][k]=new TH1D(Form("hfRing2_%d_%d_%d",i,j,k),Form("hfRing2_%d_%d_%d",i,j,k),100,0.,1.);
			   hfRing2[i][j][k]->SetTitle(";frac03;");
			   hfRing3[i][j][k]=new TH1D(Form("hfRing3_%d_%d_%d",i,j,k),Form("hfRing3_%d_%d_%d",i,j,k),100,0.,1.);
			   hfRing3[i][j][k]->SetTitle(";frac04;");
			   hptD[i][j][k]=new TH1D(Form("hptD_%d_%d_%d",i,j,k),Form("hptD_%d_%d_%d",i,j,k),100,0.,1.);
			   hptD[i][j][k]->SetTitle(";ptD;");
			   hbeta[i][j][k]=new TH1D(Form("hbeta_%d_%d_%d",i,j,k),Form("hbeta_%d_%d_%d",i,j,k),100,0.,1.);
			   hbeta[i][j][k]->SetTitle(";beta;");
			   hpull[i][j][k]=new TH1D(Form("hpull_%d_%d_%d",i,j,k),Form("hpull_%d_%d_%d",i,j,k),25,0.,0.025);
			   hpull[i][j][k]->SetTitle(";pull;");
			   hjetR[i][j][k]=new TH1D(Form("hjetR_%d_%d_%d",i,j,k),Form("hjetR_%d_%d_%d",i,j,k),100,0.,1.);
			   hjetR[i][j][k]->SetTitle(";jetR;");
			   hjetRchg[i][j][k]=new TH1D(Form("hjetRchg_%d_%d_%d",i,j,k),Form("hjetRchg_%d_%d_%d",i,j,k),100,0.,1.);
			   hjetRchg[i][j][k]->SetTitle(";jetRchg;");
			   hjtpt[i][j][k]=new TH1D(Form("hjtpt_%d_%d_%d",i,j,k),Form("hjtpt_%d_%d_%d",i,j,k),100,0.,100.);
			   hjtpt[i][j][k]->SetTitle(";jet p_{T};");
			   hjteta[i][j][k]=new TH1D(Form("hjteta_%d_%d_%d",i,j,k),Form("hjteta_%d_%d_%d",i,j,k),100,-5.,5.);
			   hjteta[i][j][k]->SetTitle(";jet #eta;");
 
		   }
	   }
   }

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"total entries: "<<nentries<<endl;

   bool pre_sel;
   bool flavor_def[3];
   bool pt_bin[3];
   bool eta_bin[4];
   float myweight;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
	  if (jentry%1000==0)
		cout<<"i = "<<jentry<<endl;
	  if(weight>0)
		myweight=1.;
	  else
		myweight=-1.;
	  for(unsigned int index=0;index<DRweighted->size();index++){
		  pre_sel=(jtpt->at(index)>20 && jtpt->at(index)<100 && abs(jteta->at(index))<5. && b_atLeastOneMuon_TrigMatch==1 && dR_Jet_Muon->at(index)>0.4);;
		  flavor_def[0]=(refdrjt->at(index)<0.2 && abs(partonFlavor->at(index))>0 && abs(partonFlavor->at(index))<6); // quarks
		  flavor_def[1]=(refdrjt->at(index)<0.2 && partonFlavor->at(index)==21); // glouns
		  flavor_def[2]=(refdrjt->at(index)>0.3 && partonFlavor->at(index)==0); // pileup
		  pt_bin[0]=(jtpt->at(index)>20 && jtpt->at(index)<30);
		  pt_bin[1]=(jtpt->at(index)>30 && jtpt->at(index)<50);
		  pt_bin[2]=(jtpt->at(index)>50 && jtpt->at(index)<100);
		  eta_bin[0]=(abs(jteta->at(index))<2.5);
		  eta_bin[1]=(abs(jteta->at(index))>2.5 && abs(jteta->at(index))<2.75);
		  eta_bin[2]=(abs(jteta->at(index))>2.75 && abs(jteta->at(index))<3.);
		  eta_bin[3]=(abs(jteta->at(index))>3. && abs(jteta->at(index))<5.);
		  if(pre_sel){
			  for(int i=0;i<3;i++){
				  if(!flavor_def[i]) continue;
				  for(int j=0;j<3;j++){
					  if(!pt_bin[j]) continue;
					  for(int k=0;k<4;k++){
						  if(!eta_bin[k]) continue;
						  hDRweighted[i][j][k]->Fill(DRweighted->at(index),myweight);
						  if(index==0)
							hnpv[i][j][k]->Fill(npv,myweight);
						  hnTot[i][j][k]->Fill(nTot->at(index),myweight);
						  hnCh[i][j][k]->Fill(nCh->at(index),myweight);
						  haxisMajor[i][j][k]->Fill(axisMajor->at(index),myweight);
						  haxisMinor[i][j][k]->Fill(axisMinor->at(index),myweight);
						  hfRing0[i][j][k]->Fill(fRing0->at(index),myweight);
						  hfRing1[i][j][k]->Fill(fRing1->at(index),myweight);
						  hfRing2[i][j][k]->Fill(fRing2->at(index),myweight);
						  hfRing3[i][j][k]->Fill(fRing3->at(index),myweight);
						  hptD[i][j][k]->Fill(ptD->at(index),myweight);
						  hbeta[i][j][k]->Fill(beta->at(index),myweight);
						  hpull[i][j][k]->Fill(pull->at(index),myweight);
						  hjetR[i][j][k]->Fill(jetR->at(index),myweight);
						  hjetRchg[i][j][k]->Fill(jetRchg->at(index),myweight);
						  hjtpt[i][j][k] ->Fill(jtpt->at(index),myweight);
						  hjteta[i][j][k]->Fill(jteta->at(index),myweight);

					  }
				  }
			  }
		  }
	  }
   }

   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(0);

   TCanvas * cDRweighted = new TCanvas("cDRweighted","cDRweighted",900,900);
   cDRweighted->Divide(3,4);
   int icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cDRweighted->cd(icanvas);
		   hDRweighted[0][j][k]->SetLineColor(kRed);
		   hDRweighted[1][j][k]->SetLineColor(kBlue);
		   hDRweighted[2][j][k]->SetLineColor(kGreen);
		   hDRweighted[0][j][k]->DrawNormalized("histSame");
		   hDRweighted[1][j][k]->DrawNormalized("histSame");
		   hDRweighted[2][j][k]->DrawNormalized("histSame");
       }
   }
   cDRweighted->Print("DRweighted.png");

   TCanvas * cnpv = new TCanvas("cnpv","cnpv",900,900);
   cnpv->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cnpv->cd(icanvas);
		   hnpv[0][j][k]->SetLineColor(kRed);
		   hnpv[1][j][k]->SetLineColor(kBlue);
		   hnpv[2][j][k]->SetLineColor(kGreen);
		   hnpv[0][j][k]->DrawNormalized("histSame");
		   hnpv[1][j][k]->DrawNormalized("histSame");
		   hnpv[2][j][k]->DrawNormalized("histSame");
       }
   }
   cnpv->Print("npv.png");

   TCanvas * cnTot = new TCanvas("cnTot","cnTot",900,900);
   cnTot->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cnTot->cd(icanvas);
		   hnTot[0][j][k]->SetLineColor(kRed);
		   hnTot[1][j][k]->SetLineColor(kBlue);
		   hnTot[2][j][k]->SetLineColor(kGreen);
		   hnTot[1][j][k]->DrawNormalized("histSame");
		   hnTot[0][j][k]->DrawNormalized("histSame");
		   hnTot[2][j][k]->DrawNormalized("histSame");
       }
   }
   cnTot->Print("nTot.png");

   TCanvas * cnCh = new TCanvas("cnCh","cnCh",900,900);
   cnCh->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cnCh->cd(icanvas);
		   hnCh[0][j][k]->SetLineColor(kRed);
		   hnCh[1][j][k]->SetLineColor(kBlue);
		   hnCh[2][j][k]->SetLineColor(kGreen);
		   hnCh[0][j][k]->DrawNormalized("histSame");
		   hnCh[1][j][k]->DrawNormalized("histSame");
		   hnCh[2][j][k]->DrawNormalized("histSame");
       }
   }
   cnCh->Print("nCh.png");

   TCanvas * caxisMajor = new TCanvas("caxisMajor","caxisMajor",900,900);
   caxisMajor->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   caxisMajor->cd(icanvas);
		   haxisMajor[0][j][k]->SetLineColor(kRed);
		   haxisMajor[1][j][k]->SetLineColor(kBlue);
		   haxisMajor[2][j][k]->SetLineColor(kGreen);
		   haxisMajor[0][j][k]->DrawNormalized("histSame");
		   haxisMajor[1][j][k]->DrawNormalized("histSame");
		   haxisMajor[2][j][k]->DrawNormalized("histSame");
       }
   }
   caxisMajor->Print("axisMajor.png");

   TCanvas * caxisMinor = new TCanvas("caxisMinor","caxisMinor",900,900);
   caxisMinor->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   caxisMinor->cd(icanvas);
		   haxisMinor[0][j][k]->SetLineColor(kRed);
		   haxisMinor[1][j][k]->SetLineColor(kBlue);
		   haxisMinor[2][j][k]->SetLineColor(kGreen);
		   haxisMinor[0][j][k]->DrawNormalized("histSame");
		   haxisMinor[1][j][k]->DrawNormalized("histSame");
		   haxisMinor[2][j][k]->DrawNormalized("histSame");
       }
   }
   caxisMinor->Print("axisMinor.png");

   TCanvas * cfRing0 = new TCanvas("cfRing0","cfRing0",900,900);
   cfRing0->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cfRing0->cd(icanvas);
		   hfRing0[0][j][k]->SetLineColor(kRed);
		   hfRing0[1][j][k]->SetLineColor(kBlue);
		   hfRing0[2][j][k]->SetLineColor(kGreen);
		   hfRing0[0][j][k]->DrawNormalized("histSame");
		   hfRing0[1][j][k]->DrawNormalized("histSame");
		   hfRing0[2][j][k]->DrawNormalized("histSame");
       }
   }
   cfRing0->Print("fRing0.png");

   TCanvas * cfRing1 = new TCanvas("cfRing1","cfRing1",900,900);
   cfRing1->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cfRing1->cd(icanvas);
		   hfRing1[0][j][k]->SetLineColor(kRed);
		   hfRing1[1][j][k]->SetLineColor(kBlue);
		   hfRing1[2][j][k]->SetLineColor(kGreen);
		   hfRing1[0][j][k]->DrawNormalized("histSame");
		   hfRing1[1][j][k]->DrawNormalized("histSame");
		   hfRing1[2][j][k]->DrawNormalized("histSame");
       }
   }
   cfRing1->Print("fRing1.png");

   TCanvas * cfRing2 = new TCanvas("cfRing2","cfRing2",900,900);
   cfRing2->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cfRing2->cd(icanvas);
		   hfRing2[0][j][k]->SetLineColor(kRed);
		   hfRing2[1][j][k]->SetLineColor(kBlue);
		   hfRing2[2][j][k]->SetLineColor(kGreen);
		   hfRing2[0][j][k]->DrawNormalized("histSame");
		   hfRing2[1][j][k]->DrawNormalized("histSame");
		   hfRing2[2][j][k]->DrawNormalized("histSame");
       }
   }
   cfRing2->Print("fRing2.png");

   TCanvas * cfRing3 = new TCanvas("cfRing3","cfRing3",900,900);
   cfRing3->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cfRing3->cd(icanvas);
		   hfRing3[0][j][k]->SetLineColor(kRed);
		   hfRing3[1][j][k]->SetLineColor(kBlue);
		   hfRing3[2][j][k]->SetLineColor(kGreen);
		   hfRing3[0][j][k]->DrawNormalized("histSame");
		   hfRing3[1][j][k]->DrawNormalized("histSame");
		   hfRing3[2][j][k]->DrawNormalized("histSame");
       }
   }
   cfRing3->Print("fRing3.png");

   TCanvas * cptD = new TCanvas("cptD","cptD",900,900);
   cptD->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cptD->cd(icanvas);
		   hptD[0][j][k]->SetLineColor(kRed);
		   hptD[1][j][k]->SetLineColor(kBlue);
		   hptD[2][j][k]->SetLineColor(kGreen);
		   hptD[1][j][k]->DrawNormalized("histSame");
		   hptD[0][j][k]->DrawNormalized("histSame");
		   hptD[2][j][k]->DrawNormalized("histSame");
       }
   }
   cptD->Print("ptD.png");

   TCanvas * cbeta = new TCanvas("cbeta","cbeta",900,900);
   cbeta->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cbeta->cd(icanvas);
		   hbeta[0][j][k]->SetLineColor(kRed);
		   hbeta[1][j][k]->SetLineColor(kBlue);
		   hbeta[2][j][k]->SetLineColor(kGreen);
		   hbeta[0][j][k]->DrawNormalized("histSame");
		   hbeta[1][j][k]->DrawNormalized("histSame");
		   hbeta[2][j][k]->DrawNormalized("histSame");
       }
   }
   cbeta->Print("beta.png");

   TCanvas * cpull = new TCanvas("cpull","cpull",900,900);
   cpull->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cpull->cd(icanvas);
		   hpull[0][j][k]->SetLineColor(kRed);
		   hpull[1][j][k]->SetLineColor(kBlue);
		   hpull[2][j][k]->SetLineColor(kGreen);
		   hpull[0][j][k]->DrawNormalized("histSame");
		   hpull[1][j][k]->DrawNormalized("histSame");
		   hpull[2][j][k]->DrawNormalized("histSame");
       }
   }
   cpull->Print("pull.png");

   TCanvas * cjetR = new TCanvas("cjetR","cjetR",900,900);
   cjetR->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cjetR->cd(icanvas);
		   hjetR[0][j][k]->SetLineColor(kRed);
		   hjetR[1][j][k]->SetLineColor(kBlue);
		   hjetR[2][j][k]->SetLineColor(kGreen);
		   hjetR[1][j][k]->DrawNormalized("histSame");
		   hjetR[0][j][k]->DrawNormalized("histSame");
		   hjetR[2][j][k]->DrawNormalized("histSame");
       }
   }
   cjetR->Print("jetR.png");

   TCanvas * cjetRchg = new TCanvas("cjetRchg","cjetRchg",900,900);
   cjetRchg->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cjetRchg->cd(icanvas);
		   hjetRchg[0][j][k]->SetLineColor(kRed);
		   hjetRchg[1][j][k]->SetLineColor(kBlue);
		   hjetRchg[2][j][k]->SetLineColor(kGreen);
		   hjetRchg[1][j][k]->DrawNormalized("histSame");
		   hjetRchg[0][j][k]->DrawNormalized("histSame");
		   hjetRchg[2][j][k]->DrawNormalized("histSame");
       }
   }
   cjetRchg->Print("jetRchg.png");

   TCanvas * cjtpt = new TCanvas("cjtpt","cjtpt",900,900);
   cjtpt->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cjtpt->cd(icanvas);
		   hjtpt[0][j][k]->SetLineColor(kRed);
		   hjtpt[1][j][k]->SetLineColor(kBlue);
		   hjtpt[2][j][k]->SetLineColor(kGreen);
		   hjtpt[0][j][k]->DrawNormalized("histSame");
		   hjtpt[1][j][k]->DrawNormalized("histSame");
		   hjtpt[2][j][k]->DrawNormalized("histSame");
       }
   }
   cjtpt->Print("jtpt.png");

   TCanvas * cjteta = new TCanvas("cjteta","cjteta",900,900);
   cjteta->Divide(3,4);
   icanvas=0;
   for(int k=0;k<4;k++){
       for(int j=0;j<3;j++){
		   icanvas++;
		   cjteta->cd(icanvas);
		   hjteta[0][j][k]->SetLineColor(kRed);
		   hjteta[1][j][k]->SetLineColor(kBlue);
		   hjteta[2][j][k]->SetLineColor(kGreen);
		   hjteta[0][j][k]->DrawNormalized("histSame");
		   hjteta[1][j][k]->DrawNormalized("histSame");
		   hjteta[2][j][k]->DrawNormalized("histSame");
       }
   }
   cjteta->Print("jteta.png");

}

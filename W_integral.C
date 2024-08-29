#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include <TStopwatch.h>

//Delcare Global Parameters
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event

void W_integral(Int_t entries_input = -1,Int_t kine = 2) {// Main

 TChain *C = new TChain("T");

 Double_t HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
 Double_t HCal_th = 35.0; // Angle that the center of HCal is at 

 string configfilename = Form("/w/halla-scshelf2102/sbs/cmjackso/GEn/config/W2_config.cfg");
  cout<<"reading in from a config file"<<endl;

  ifstream configfile(configfilename);
 TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ) {
    if( !currentline.BeginsWith("#") ) {
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      C->Add(currentline);
      cout << "Loaded file at: " << currentline << endl;
    }    
  }
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ) {
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }    
    cout<< "Global Cut: "<<globalcut<<endl;
  }
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ) {
    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
	cout << "Loading HCal distance: " << HCal_d << endl;
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof();	
	cout << "Loading HCal angle: " << HCal_th << endl;
      }
    }
    delete tokens;
  }

  //Apply Global Cut
  cout<<endl<<"Populating list with global cut. May take a few minutes...."<<endl;
  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);
  cout << endl << "Event list populated with cut: "<<globalcut << endl;

   // Set long int to keep track of total entries
  Long64_t Nevents = elist->GetN();
  UInt_t run_number = 0;

  cout<<endl << "Opened up TChain with nentries: " << C->GetEntries() << ". After globalcut: " << Nevents << "." << endl << endl;
  cout<<"Entries: "<<Nevents<<endl;
 

 //Declaring Parameters         
 // Double_t shIdblk;
 Double_t runnum[2100];
 // HCAL
 Double_t HCAL_e;
 // BBCAL
  Double_t BBps_e, BBsh_e;
  Double_t bb_tdctrig_tdc[1000], bb_tdctrig_tdcelemID[1000];
  Int_t Nclusters; 
  Int_t Ndata_bb_tdctrig_tdcelemID;
 // Track params 
  Double_t BBtr_x[maxTracks], BBtr_y[maxTracks], BBtr_vz[maxTracks], BBtr_p[maxTracks], BBtr_th[maxTracks], BBtr_ph[maxTracks];
  Double_t BBtr_n;
  // Hodoscope
  Double_t HODOtmean[1000];
  // Physics
  Double_t W2;

  //GRINCH PARAMETERS: Turn off if not needed

   Double_t BBgr_hit_amp[1000],BBgr_hit_pmtnum[1000],BBgr_hit_time[1000],BBgr_hit_clustindex[1000];
  Int_t hitsGR;
 Double_t BBgr_allclus_tmean[1000], BBgr_allclus_adc[1000], BBgr_allclus_size[1000], BBgr_allclus_trms[1000], BBgr_allclus_tot_mean[1000], BBgr_allclus_trackindex[1000], BBgr_allclus_xmean[1000], BBgr_allclus_ymean[1000], BBgr_allclus_dx[1000], BBgr_allclus_dy[1000]; 

 Double_t BBgr_allclus_mirrorindex[1000];
  Double_t BBgr_clus_mirrorindex, BBgr_clus_tmean, BBgr_clus_adc, BBgr_clus_size, BBgr_clus_trms, BBgr_clus_tot_mean, BBgr_clus_trackindex, BBgr_clus_xmean, BBgr_clus_ymean; 
  Double_t BBgr_bestcluster, BBgr_ngoodhits, BBgr_ntrackmatch;


  // Int_t Ndata_BBgr_tdcelemID;
  // Double_t  BBgr_tdc_tdcelemID[1000], BBgr_tdc_tdc[1000],BBgr_tdc_te[1000], BBgr_tdc_mult[1000], BBgr_tdc_tot[1000];

 // Declare root tree variabls and set values to memeory locations in root file 

 // Turn off all branches
  C->SetBranchStatus( "*", 0 ); 
 // Turn on only branches we need.
  C->SetBranchStatus( "g.runnum", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.x", 1 );
  C->SetBranchStatus( "bb.tr.y", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "bb.tr.vz", 1 );  
  C->SetBranchStatus( "bb.tr.th", 1 );
  C->SetBranchStatus( "bb.tr.ph", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
  C->SetBranchStatus("bb.tdctrig.tdc",1);
  C->SetBranchStatus("bb.tdctrig.tdcelemID",1);
  C->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);
  C->SetBranchStatus("e.kine.W2",1);

  //GRINCH VARIABLES: Turn off when not needed

   //All-clust variables
  C->SetBranchStatus("Ndata.bb.grinch_tdc.allclus.adc",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.adc",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.size",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.t_mean",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.t_rms",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.tot_mean",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.trackindex",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.x_mean",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.y_mean",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.mirrorindex",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.dx",1);//
  C->SetBranchStatus("bb.grinch_tdc.allclus.dy",1);//

  //track matched cluster
  C->SetBranchStatus("bb.grinch_tdc.clus.adc",1); // TDC LE sum 
  C->SetBranchStatus("bb.grinch_tdc.clus.size",1); // number of PMTs in the cluster
  C->SetBranchStatus("bb.grinch_tdc.clus.t_mean",1); // average LE time pf the PMTs
  C->SetBranchStatus("bb.grinch_tdc.clus.t_rms",1); // RMS of the average LE time. 
  C->SetBranchStatus("bb.grinch_tdc.clus.tot_mean",1); // mean of the time-over-threshold
  C->SetBranchStatus("bb.grinch_tdc.clus.trackindex",1);// which track the cluster matches (-1 if none)
  C->SetBranchStatus("bb.grinch_tdc.clus.x_mean",1); // the mean x position of the PMTs in the cluster
  C->SetBranchStatus("bb.grinch_tdc.clus.y_mean",1); //  the mean y position of the PMTs in the cluster
  C->SetBranchStatus("bb.grinch_tdc.clus.mirrorindex",1); // which mirror it was matched to
  C->SetBranchStatus("bb.grinch_tdc.hit.clustindex",1);
  C->SetBranchStatus("bb.grinch_tdc.bestcluster",1);
  C->SetBranchStatus("bb.grinch_tdc.ngoodhits",1);
  C->SetBranchStatus("bb.grinch_tdc.ntrackmatch",1);
  C->SetBranchStatus( "bb.grinch_tdc.hit.amp", 1 );
  C->SetBranchStatus( "bb.grinch_tdc.hit.time", 1 );
  C->SetBranchStatus( "bb.grinch_tdc.hit.pmtnum", 1 );
  C->SetBranchStatus( "Ndata.bb.grinch_tdc.hit.pmtnum", 1 );
  

   // Map branches to the variables 
  C->SetBranchAddress( "sbs.hcal.e", &HCAL_e );
  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
  C->SetBranchAddress( "bb.tr.x", &BBtr_x );
  C->SetBranchAddress( "bb.tr.y", &BBtr_y );
  C->SetBranchAddress( "bb.tr.vz", &BBtr_vz ); 
  C->SetBranchAddress( "bb.tr.p", &BBtr_p );
  C->SetBranchAddress( "bb.tr.th", &BBtr_th );
  C->SetBranchAddress( "bb.tr.ph", &BBtr_ph );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean );
  C->SetBranchAddress("bb.tdctrig.tdc",&bb_tdctrig_tdc);
  C->SetBranchAddress("bb.tdctrig.tdcelemID",&bb_tdctrig_tdcelemID);
  C->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID",&Ndata_bb_tdctrig_tdcelemID);
  C->SetBranchAddress( "e.kine.W2", &W2 );
  C->SetBranchAddress( "g.runnum",&runnum); // Run Number

  //GRINCH MAPPED VARIABLES: Turn off when not needed
  C->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.adc",&Nclusters);
  C->SetBranchAddress("bb.grinch_tdc.allclus.adc",&BBgr_allclus_adc);
  C->SetBranchAddress("bb.grinch_tdc.allclus.size",&BBgr_allclus_size);
  C->SetBranchAddress("bb.grinch_tdc.allclus.t_mean",&BBgr_allclus_tmean);
  C->SetBranchAddress("bb.grinch_tdc.allclus.t_rms",&BBgr_allclus_trms);
  C->SetBranchAddress("bb.grinch_tdc.allclus.tot_mean",&BBgr_allclus_tot_mean);
  C->SetBranchAddress("bb.grinch_tdc.allclus.trackindex",&BBgr_allclus_trackindex);
  C->SetBranchAddress("bb.grinch_tdc.allclus.x_mean",&BBgr_allclus_xmean);
  C->SetBranchAddress("bb.grinch_tdc.allclus.y_mean",&BBgr_allclus_ymean);
  C->SetBranchAddress("bb.grinch_tdc.allclus.mirrorindex",&BBgr_allclus_mirrorindex);
  C->SetBranchAddress("bb.grinch_tdc.allclus.dx",&BBgr_allclus_dx);//
  C->SetBranchAddress("bb.grinch_tdc.allclus.dy",&BBgr_allclus_dy);//

  C->SetBranchAddress( "bb.grinch_tdc.hit.amp",&BBgr_hit_amp); // ToT
  C->SetBranchAddress( "bb.grinch_tdc.hit.time",&BBgr_hit_time); //Leading Edge
  C->SetBranchAddress( "bb.grinch_tdc.hit.pmtnum",&BBgr_hit_pmtnum );
  C->SetBranchAddress( "Ndata.bb.grinch_tdc.hit.pmtnum",&hitsGR ); // Num PMT fired

  C->SetBranchAddress("bb.grinch_tdc.hit.clustindex",&BBgr_hit_clustindex);
  C->SetBranchAddress("bb.grinch_tdc.clus.adc",&BBgr_clus_adc);
  C->SetBranchAddress("bb.grinch_tdc.clus.size",&BBgr_clus_size);
  C->SetBranchAddress("bb.grinch_tdc.clus.t_mean",&BBgr_clus_tmean);
  C->SetBranchAddress("bb.grinch_tdc.clus.t_rms",&BBgr_clus_trms);
  C->SetBranchAddress("bb.grinch_tdc.clus.tot_mean",&BBgr_clus_tot_mean);
  C->SetBranchAddress("bb.grinch_tdc.clus.trackindex",&BBgr_clus_trackindex);
  C->SetBranchAddress("bb.grinch_tdc.clus.x_mean",&BBgr_clus_xmean);
  C->SetBranchAddress("bb.grinch_tdc.clus.y_mean",&BBgr_clus_ymean);
  C->SetBranchAddress("bb.grinch_tdc.clus.mirrorindex",&BBgr_clus_mirrorindex);
  C->SetBranchAddress("bb.grinch_tdc.bestcluster",&BBgr_bestcluster);
  C->SetBranchAddress("bb.grinch_tdc.ngoodhits",&BBgr_ngoodhits);
  C->SetBranchAddress("bb.grinch_tdc.ntrackmatch",&BBgr_ntrackmatch);

  
  //Declare OutFiles
 TFile *fout = new TFile( Form("/w/halla-scshelf2102/sbs/cmjackso/GEN_analysis/basic/outfiles/W2_looky_trackcuts_REPLAY_again.root"), "RECREATE" );

 //Create New TTree
 TTree *outTree =C->CloneTree(0);
 
 //Define histograms

 TH1D* h_W2 = new TH1D("h_W2",";W2",200,0.0,3.0);
 //h_W2->GetXaxis()->SetTitle("W^2 [GeV^2]");

 //TH1D* h_W_loss = new TH1D("h_W_loss",";W without Low Signal Shower Channels",150,0.0,3.0);
 //h_W_loss->GetXaxis()->SetTitle("W^2 [GeV^2]");
 TH1D* h_BBps_e =  new TH1D("h_BBps_e",";BBps_e;" ,200, 0, 2);
 TH1D* h_HCAL_e = new TH1D("h_HCAL_e",";HCAL e;",200,0,3);
 TH1D* h_BBsh_e =  new TH1D("h_BBsh_e",";BBsh_e;" ,200, 0, 3);
 TH2D* h_BBsh_ps_e = new TH2D("h_BBsh_ps_e", "; sh_e  ; ps_e ", 100,0,3,100,0,3 );
 TH1D* h_BB_e_p =  new TH1D("h_BB_e_p",";total shower energy/ p;" ,200, 0, 2);


 
  
 //check if input for the number of events is valid
  Int_t max = 0;
  if (entries_input == -1){
    max = Nevents;
  }
  else{
    max = entries_input;
  }
  if (max > Nevents){ max = Nevents;}
  cout<<"max = "<<max<<endl; 

  // Loop over events
  for(Long64_t nevent = 0; nevent<max; nevent++) {
   
    C->GetEntry( elist->GetEntry( nevent ) ); 
    if (nevent%100000==0) cout << " Entry = " << nevent << endl;

    Double_t projx = BBtr_x[0]+BBtr_th[0]*0.48;//0.48
    Double_t projy = BBtr_y[0]+BBtr_ph[0]*0.48;

    // General Histos
    h_W2 ->Fill(W2);
    h_HCAL_e ->Fill(HCAL_e);
    h_BBps_e ->Fill(BBps_e);
    h_BBsh_e ->Fill(BBsh_e);
    h_BBsh_ps_e ->Fill(BBsh_e,BBps_e);
    Double_t e_over_p  = (BBps_e + BBsh_e)/BBtr_p[0];
    h_BB_e_p -> Fill(e_over_p);

    //Fill the New TTree
    outTree->Fill();
  }//end loop over entries
  //Write histograms & TTree to File
    fout->cd();
    outTree->Write();
    h_W2->Write();
    h_BBps_e->Write();
    h_HCAL_e->Write();
    h_BBsh_e->Write();
    h_BBsh_ps_e->Write();
    h_BB_e_p->Write();
    fout->Close();
 // TCanvas *c1 = new TCanvas("c1","Test",100,100,700,700);
 // c1->Divide(1,2,3,4);
 // c1->cd(1);
 // h_W2->Draw();
 // c1->cd(2);
 // h_BBps_e->Draw();
    // c1->cd(3);
    
    cout << "Data merged and written to " << fout->GetName() << endl;
    //fout -> Write();
 
}
  //end main

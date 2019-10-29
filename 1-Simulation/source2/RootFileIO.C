#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TKey.h"
#include "iostream"

//*****************************************************************************
// To run this macro in root:
//
// gSystem->Load("/Users/liuqian/Work/work/geant4/Geant4 Framework/v2/debug/libSimEventDict.so");
// .L RootFileIO.C++
// RootFileIO();
//
//*****************************************************************************

#include "/home/stuf/geant4.10.03.p03/examples/basic/geant4track/source/include/MyHistMathLibrary.h"
#include "/home/stuf/geant4.10.03.p03/examples/basic/geant4track/source/include/SimEvent.h"

void WriteRoot();
void ReadRoot();

TString fFileName("tmp.root");

void RootFileIO()
{

//   WriteRoot();

   ReadRoot();

}


void WriteRoot()
{
  cout<<"---> Write to "<<fFileName<<" root file."<<endl;

  TFile * fRootFp = new TFile(fFileName, "recreate");

  SimEvent* fEvent = new SimEvent();
  TTree   * fTree1 = new TTree("sim", "Tree of GasTPC events");


  fTree1->Branch("SimEvent", "SimEvent", &fEvent, 32000, 99);
  //for(int i=0; i<9; i++) fEvent->SetTime(i, 10*i);

  fTree1->Fill();
  fEvent->Clear();

  fRootFp->Write();
  fRootFp->Close();
}

void ReadRoot()
{
  TH2F *f = new TH2F("f", "f", 10, 1, 10, 10, 1, 10);
  for(int i=0; i<100; i++) f->Fill(5, 5);
  
  TH2F *f2 = new TH2F("f2", "f", 10, 1, 10, 10, 1, 10);
  SmearingGaus(f, f2, 0.5);

  return;
  cout<<"---> Read from "<<fFileName<<" root file."<<endl;

  TFile * fRootFp = new TFile(fFileName);    
  
  SimEvent* fEvent = 0;
  TTree   * fTree1  = (TTree *)fRootFp->Get("sim");
  fTree1->SetBranchAddress("SimEvent", &fEvent);

  fTree1->GetEntry(0);
  //cout<<fEvent->GetTime(2)<<endl;  

  fRootFp->Close();
}



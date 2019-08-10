// ********************************************************************
//
#include "TApplication.h"
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TKey.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TView.h"
#include "TPolyLine3D.h"
//
#include "SimEvent.h"
#include "G4ios.hh"
#include <iostream>
#include "Riostream.h"
#include "MyHistMathLibrary.h"

#include "time.h"

void WriteRoot(TString fFileName);
void ReadRoot(TString fFileName);

int verbose = 2;

int main(int argc, char **argv)
{
    TApplication theApp("App", &argc, argv); 

    TString fFileName = TString("tmp.root");

    // initialize ROOT
    TSystem ts;
    gSystem->Load("libSimEventDict");

    if (argc < 2)
        G4cout << "Missing name of the file to read! Read " << fFileName << " instead." << G4endl;
    else
        fFileName = argv[1]+TString(".root");  //自动添加后缀名

    G4cout << " Read " << fFileName << G4endl;
    ReadRoot(fFileName);
    theApp.Run();  


    return 0;
}

void WriteRoot(TString fFileName)
{
    if (verbose >= 0)
        G4cout << "---> Write to " << fFileName << " root file." << G4endl;

    TFile *fRootFp = new TFile(fFileName, "recreate");

    SimEvent *fEvent = new SimEvent();
    TTree *fTree1 = new TTree("sim", "Tree of GasTPC events");

    fTree1->Branch("SimEvent", "SimEvent", &fEvent, 32000, 99);

    fTree1->Fill();
    fEvent->Clear();

    fRootFp->Write();
    fRootFp->Close();
}



void ReadRoot(TString fFileName)
{
    if (verbose >= 0)
        G4cout << "---> Reading from " << fFileName << " root file." << G4endl;

    //---- draw histogram ----
    //
    const int NX = 72;
    const int NY = 72;
    const double lenPit = 0.08;
    const double conFac = 10000 / 30E-6; //gain/W_avg*(??DAQ??)

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 480);
    TH2F *fhitXY = new TH2F("fhitXY", "fhitXY", NX, -lenPit * 36, lenPit * 36, NY, -lenPit * 36, lenPit * 36);
    TH2F *fhitXZ = new TH2F("fhitXZ", "fhitXZ", NX, -lenPit * 36, lenPit * 36, NY, -lenPit * 72, 0);
    TH2F *fhitYZ = new TH2F("fhitYZ", "fhitYZ", NX, -lenPit * 36, lenPit * 36, NY, -lenPit * 72, 0);

    //---- init root file to read ----
    //
    TFile *fRootFp = new TFile(fFileName);
    if (!fRootFp->IsOpen())
    {
        G4cout << "---> Can't open " << fFileName << "." << G4endl;
        return;
    }

    SimEvent *fEvent = 0;
    TTree *fTree1 = (TTree *)fRootFp->Get("sim");
    fTree1->SetBranchAddress("SimEvent", &fEvent);

    //---- process data-reading ----
    //
    double rates = fTree1->GetEntries();
    G4cout << "> There are " << rates << " events stored." << G4endl;

    for (int i = 0; i < 1; i++)//fTree1->GetEntries(); i++)
    {
        if (verbose >= 0 && i % 100 == 0)
            G4cout << "> No. " << i << " event: " << G4endl;

        // 0. read entry	
        fTree1->GetEntry(i);
        fhitXY->Reset();
        fhitXZ->Reset();
        fhitYZ->Reset();

        // 1. the global information of the fEvent
        if (verbose > 0)
        {
            G4cout << ">  PdgID : " << fEvent->GetPDGID() << G4endl;
            G4cout << ">  Energy: " << fEvent->GetTrueEnergy() << G4endl;
            G4cout << ">  Moment: (" << fEvent->GetMomentumDir().x() << ", "
                   << fEvent->GetMomentumDir().y() << ", "
                   << fEvent->GetMomentumDir().z() << ") " << G4endl;
        }

        // 2. read tracks information
        TMap *trkMap1 = fEvent->GetGasSDTrackMap();
        TMap *depMap1 = fEvent->GetGasSDDepositMap();

        TMap *trkMap2 = fEvent->GetWindowSDTrackMap();
        TMap *depMap2 = fEvent->GetWindowSDDepositMap();

        if (verbose)
            G4cout << ">-----Tracks in Gas : " << fEvent->GetGasDepEnergy() << " MeV.\n";

        TObject *tempObj = 0;
        TIterator *it = trkMap1->MakeIterator();
        while ((tempObj = it->Next()))
        {
            TObjString *obj = dynamic_cast<TObjString *>(tempObj);
            SimTrack *trk = dynamic_cast<SimTrack *>(trkMap1->GetValue(obj));
            if (verbose > 1)
            {
                G4cout << ": \t-------------------------" << G4endl;
                G4cout << ": \t-PdgID :" << trk->GetPDGID() << G4endl;
                G4cout << ": \t TrkID :" << trk->GetTrackID() << G4endl;
                G4cout << ": \t Mass  :" << trk->GetInitMass() << G4endl;
                G4cout << ": \t Ek    :" << trk->GetInitEk() << G4endl;
                G4cout << ": \t Edep  :" << trk->GetEdep() << G4endl;
                G4cout << ": \t IniPx :" << trk->GetInitPx() << G4endl;
                G4cout << ": \t IniPy :" << trk->GetInitPy() << G4endl;
                G4cout << ": \t IniPz :" << trk->GetInitPz() << G4endl;
                G4cout << ": \t IniX  :" << trk->GetInitX() << G4endl;
                G4cout << ": \t IniY  :" << trk->GetInitY() << G4endl;
                G4cout << ": \t IniZ  :" << trk->GetInitZ() << G4endl;
                G4cout << ": \t IniT  :" << trk->GetInitT() << G4endl;
                G4cout << ": \t ExtPx :" << trk->GetExitPx() << G4endl;
                G4cout << ": \t ExtPy :" << trk->GetExitPy() << G4endl;
                G4cout << ": \t ExtPz :" << trk->GetExitPz() << G4endl;
                G4cout << ": \t ExtX  :" << trk->GetExitX() << G4endl;
                G4cout << ": \t ExtY  :" << trk->GetExitY() << G4endl;
                G4cout << ": \t ExtZ  :" << trk->GetExitZ() << G4endl;
                G4cout << ": \t ExtT  :" << trk->GetExitT() << G4endl;
                G4cout << ": \t Length:" << trk->GetTrackLength() << G4endl;
            }

            std::vector<Int_t> stepIdx = trk->GetStepIdx();
            if (verbose) 
                G4cout << ": \t-Deposits in track : Particle =" << trk->GetPDGID()
                       << ", TrackID = " << trk->GetTrackID()
                       << ", ParentID = " << trk->GetParentID()
                       << ", NHits = " << stepIdx.size() << "." << G4endl;
            

            for (int j = 0; j < (int)stepIdx.size(); j++)
            {
                SimDeposit *dep = dynamic_cast<SimDeposit *>(depMap1->GetValue(Form("%d", stepIdx[j])));
                fhitXY->Fill(dep->GetPostX(), dep->GetPostY(), dep->GetEdep() * conFac);
                fhitXZ->Fill(dep->GetPostX(), dep->GetPostZ(), dep->GetEdep() * conFac);
                fhitYZ->Fill(dep->GetPostY(), dep->GetPostZ(), dep->GetEdep() * conFac);

                if (verbose > 1)
                {
                    G4cout << ": \t\t----hit  : " << j << G4endl;
                    G4cout << ": \t\t  PdgID  : " << dep->GetPDGID() << G4endl;
                    G4cout << ": \t\t  TrkID  : " << dep->GetTrackID() << G4endl;
                    G4cout << ": \t\t  ProdID : " << dep->GetProducerID() << G4endl;
                    G4cout << ": \t\t  IniPx  : " << dep->GetPrePx() << G4endl;
                    G4cout << ": \t\t  IniPy  : " << dep->GetPrePy() << G4endl;
                    G4cout << ": \t\t  IniPz  : " << dep->GetPrePz() << G4endl;
                    G4cout << ": \t\t  IniX   : " << dep->GetPreX() << G4endl;
                    G4cout << ": \t\t  IniY   : " << dep->GetPreY() << G4endl;
                    G4cout << ": \t\t  IniZ   : " << dep->GetPreZ() << G4endl;
                    G4cout << ": \t\t  IniT   : " << dep->GetPreT() << G4endl;
                    G4cout << ": \t\t  ExtPx  : " << dep->GetPostPx() << G4endl;
                    G4cout << ": \t\t  ExtPy  : " << dep->GetPostPy() << G4endl;
                    G4cout << ": \t\t  ExtPz  : " << dep->GetPostPz() << G4endl;
                    G4cout << ": \t\t  ExtX   : " << dep->GetPostX() << G4endl;
                    G4cout << ": \t\t  ExtY   : " << dep->GetPostY() << G4endl;
                    G4cout << ": \t\t  ExtZ   : " << dep->GetPostZ() << G4endl;
                    G4cout << ": \t\t  ExtT   : " << dep->GetPostT() << G4endl;
                    G4cout << ": \t\t  Edep   : " << dep->GetEdep() << G4endl;
                    G4cout << ": \t\t  Length : " << dep->GetStepLength() << G4endl;
                }
            }
        }
    }

    fRootFp->Close();

    c1->Divide(2, 2);
    c1->cd(1);
    fhitXY->Draw("colz");
    c1->cd(2);
    fhitYZ->Draw("colz");
    c1->cd(3);
    fhitXZ->Draw("colz");
    c1->Modified();
    c1->Update();
}

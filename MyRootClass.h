#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"

#include "MyEventClass.h"

const int NX = 72;
const int NY = 72;

//-----------------------------------------------------------------------------------------------
//
class MyRootClass
{
  public:
    MyRootClass(TString fp, TString pp)
    {
        ip = 0;
        nped = 3;
        hPed = NULL;
        hPedm = NULL;
        hPeds = NULL;
        hAll1 = NULL;
        filePath = fp;
        pedPath = pp;
    };
    ~MyRootClass();

    //1---
    void Analysis(TString);
    void DrawRaw(const char *opt = "colz");
    void DrawPre(const char *opt = "colz");
    void DrawNext(const char *opt = "colz");
    void DrawSelected(int, const char *opt = "colz");
    void ButtonFunc11();
    TString *GetInfo();

    //2---
    void ButtonFunc21(int, const char *opt = "colz");
    void ButtonFunc22(int, const char *opt = "colz");
    void ButtonFunc23(int, const char *opt = "colz");

    //3---
    void ButtonFunc31(double);
    void ButtonFunc32();
    void ButtonFunc33();
    void AnalysisPed(TString path);

    int GetIp() {return ip;}
    int GetNEvent() {return fEventList.size();}

  private:
    int ip;
    double nped;
    TH2F *hPed;
    TH2F *hPedm;
    TH2F *hPeds;

    TString filePath;
    TString rootFilePath;
    TString pedPath;
    TString rootPedPath;

    TH2F *hAll1; // plot all hits in one histogram
    TH2F *hAll2; // plot all hits in one histogram (remove ped)
    TH2F *hBaryCenter; // plot the all barycenter
    TH2F *hIPoint;     // plot the all IP point;



    vector<TH1F *> fPed;
    vector<MyEventClass *> fEventList;
};

MyRootClass::~MyRootClass()
{
    // Destructor.
}

//______________________________________________________________________________
void MyRootClass::AnalysisPed(TString pp)
{
    pedPath = pp;

    if (hPedm != NULL)
    {
        delete hPed;
        delete hPedm;
        delete hPeds;
    }

    hPed = new TH2F("hPed", "Pedestal", NX, 1, NX + 1, NY, 1, NY + 1);
    hPedm = new TH2F("hPedm", "Pedestal mean", NX, 1, NX + 1, NY, 1, NY + 1);
    hPeds = new TH2F("hPeds", "Pedestal sigma", NX, 1, NX + 1, NY, 1, NY + 1);

    //----------
    cout << "--> Opening: " << pedPath << endl;
    vector<int> dvec[NX][NY];
    ifstream ifSignal(pedPath, ios::binary);

    while (ifSignal.good())
    {
        unsigned short _data[NX][NY];
        ifSignal.read((char *)(&_data), sizeof(_data));
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                dvec[i][j].push_back(_data[i][j]);
    }

    //----------
    fPed.clear();
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            int xmin = dvec[i][j][0];
            int xmax = dvec[i][j][0];
            for (int k = 0; k < dvec[i][j].size(); k++)
            {
                xmin = (xmin < dvec[i][j][k]) ? xmin : dvec[i][j][k];
                xmax = (xmax > dvec[i][j][k]) ? xmax : dvec[i][j][k];
            }

            if (gDirectory->Get(Form("h%d", i * NY + j)) != NULL)
            {
                TH1F *f = (TH1F *)gDirectory->Get(Form("h%d", i * NY + j));
                delete f;
            }
            TH1F *f = new TH1F(Form("h%d", i * NY + j), Form("histogram for %d,%d", i, j), xmax - xmin + 10, xmin - 5, xmax + 5);
            for (int k = 0; k < dvec[i][j].size(); k++)
                f->Fill(dvec[i][j][k]);

            fPed.push_back(f);
        }
    }

    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            hPedm->Fill(i + 1, j + 1, fPed[i * NY + j]->GetMean());
            hPeds->Fill(i + 1, j + 1, fPed[i * NY + j]->GetRMS());
            hPed->SetBinContent(i + 1, j + 1, hPedm->GetBinContent(i + 1, j + 1) + nped * hPeds->GetBinContent(i + 1, j + 1));
        }
    }

    ifSignal.close();

    rootPedPath = pedPath;
    if (rootPedPath.Index(".mdat") != -1)
        rootPedPath.Replace(rootPedPath.Index(".mdat"), 5, ".root");
    if (rootPedPath.Index(".dat") != -1)
        rootPedPath.Replace(rootPedPath.Index(".dat"), 4, ".root");

    TFile *rootPedFile = new TFile(rootPedPath, "recreate");
    hPedm->Write();
    hPeds->Write();
    hPed->Write();
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            fPed[i * NY + j]->Write();
        }
    }
    rootPedFile->Close();

    cout << "---->Save ped root file to: " << rootPedPath << endl;

    hPedm->Draw("colz");
}

void MyRootClass::ButtonFunc31(double _ped)
{
    if (hPedm == NULL)
        return;
    if (hPeds == NULL)
        return;
    if (hPed == NULL)
        return;

    nped = _ped;

    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            hPed->SetBinContent(i + 1, j + 1, hPedm->GetBinContent(i + 1, j + 1) + nped * hPeds->GetBinContent(i + 1, j + 1));

    hPed->Draw("colz");
}

void MyRootClass::ButtonFunc32()
{
    if (hPedm != NULL)
        hPedm->Draw("colz");
}

void MyRootClass::ButtonFunc33()
{
    if (hPeds != NULL)
        hPeds->Draw("colz");
}

//______________________________________________________________________________
//
void MyRootClass::Analysis(TString filePath)
{
    if (hPedm == NULL)
    {
        rootPedPath = pedPath;
        if (rootPedPath.Index(".mdat") != -1)
            rootPedPath.Replace(rootPedPath.Index(".mdat"), 5, ".root");
        if (rootPedPath.Index(".dat") != -1)
            rootPedPath.Replace(rootPedPath.Index(".dat"), 4, ".root");
        cout << "----> Read ped from " << rootPedPath << endl;

        TFile *rootPedFile = new TFile(rootPedPath);
        hPedm = (TH2F *)rootPedFile->Get("hPedm");
        hPeds = (TH2F *)rootPedFile->Get("hPeds");
        hPed = (TH2F *)rootPedFile->Get("hPeds");

        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                hPed->SetBinContent(i + 1, j + 1, hPedm->GetBinContent(i + 1, j + 1) + nped * hPeds->GetBinContent(i + 1, j + 1));
    }

    if (hAll1 != NULL)
    {
        delete hAll1;
        delete hAll1;
        delete hBaryCenter;
        delete hIPoint;
    }
    hAll1 = new TH2F("hAll1", "Histogram of all hits", NX, 1, NX + 1, NY, 1, NY + 1);
    hAll2 = new TH2F("hAll2", "Histogram of all hits", NX, 1, NX + 1, NY, 1, NY + 1);
    hBaryCenter = new TH2F("hBaryCenter", "Histogram of all barycenters", NX, 1, NX + 1, NY, 1, NY + 1);
    hIPoint = new TH2F("hIPoint", "Histogram of all IP points", NX, 1, NX + 1, NY, 1, NY + 1);

    fEventList.clear();

    cout << "--> Opening: " << filePath << endl;
    ifstream ifSignal(filePath, ios::binary);

    ip = 1;
    int nEvent = 0;
    while (ifSignal.good())
    {
        unsigned short _data[NX][NY];
        ifSignal.read((char *)(&_data), sizeof(_data));

        if(_data[0][0] == 0xFFFF && _data[1][1] == 0xFFFF)
        {
            ifSignal.read((char *)(&_data), sizeof(_data));
            continue;
        }

        MyEventClass *fEvent = new MyEventClass(nEvent, 0, NX-1, 0, NY-1);

        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
            {
                fEvent->Fill(i, j, _data[i][j]);
                double q = _data[i][j] - hPed->GetBinContent(i + 1, j + 1);
                hAll1->Fill(i + 1, j + 1, (q < 0) ? 0 : q);
            }

        fEvent->GenerateHist(hPed);
        fEvent->Fill2DPlot(hAll2);
        fEvent->FillBaryCenter(hBaryCenter);
        fEvent->FillIPpoint(hIPoint);
        fEventList.push_back(fEvent);

        nEvent++;
        if (nEvent % 100 == 0)
            cout << "Dealed: " << nEvent << endl;
        if (nEvent > 1000)
            break;
    }

    cout << "----> " << nEvent << " events recorded." << endl;
    ifSignal.close();
}

void MyRootClass::DrawRaw(const char *opt)
{
    fEventList.at(ip)->Get2DRawPlot()->Draw(opt);
}

void MyRootClass::DrawPre(const char *opt)
{
    if (ip <= 0)
        return;

    int _ip = ip;
    while(_ip>=1)
    {
        _ip--;
        if(fEventList.at(_ip)->GetDataQuality()==1) break;
    }
    ip = _ip;
    fEventList.at(ip)->Draw2DResult(opt);
}

void MyRootClass::DrawNext(const char *opt)
{
    if ((unsigned long)ip == fEventList.size() - 1)
        return;

    int _ip = ip;
    while(_ip < fEventList.size()-1)
    {
        _ip++;
        if(fEventList.at(_ip)->GetDataQuality()==1) break;
    }
    ip = _ip;    
    fEventList.at(ip)->Draw2DResult(opt);
}

void MyRootClass::DrawSelected(int _ip, const char *opt)
{
    if(_ip < 0 || _ip > fEventList.size() - 1) 
        return;

    ip = _ip;
    fEventList.at(ip)->Draw2DResult(opt);
}

TString *MyRootClass::GetInfo()
{
    return fEventList.at(ip)->GetInfo();
}

void MyRootClass::ButtonFunc11()
{
    fEventList.at(ip)->DrawEllipse();
}

//______________________________________________________________________________

void MyRootClass::ButtonFunc21(int flag, const char *opt)
{
    if(flag==1) hAll1->Draw(opt);
    else        hAll2->Draw(opt);
}

void MyRootClass::ButtonFunc22(int flag, const char *opt)
{
    if(flag==1) hBaryCenter->Draw(opt);
    else        hIPoint->Draw(opt);
}

void MyRootClass::ButtonFunc23(int flag, const char *opt)
{
    if(flag==1) hAll1->Draw(opt);
    else        hAll2->Draw(opt);
}

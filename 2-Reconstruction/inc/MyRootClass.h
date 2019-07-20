#ifndef _MyRootClass_h_
#define _MyRootClass_h_

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGText.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TEnv.h"
#include "TSystem.h"
#include "Riostream.h"

#include "MyEventClass.h"

const int NX = 72;
const int NY = 72;

using namespace std;

//-----------------------------------------------------------------------------------------------
//
class MyRootClass
{
public:
    MyRootClass(TString fp, TString pp);
    ~MyRootClass();

    void ReadPed();
    void IniHistogram();

    //1---
    void Analysis(TString);
    void AnalysisDir(TString);
    void AnalysisFile(TString);
    void DrawRaw(const char *opt = "colz");
    void DrawPre(const char *opt = "colz");
    void DrawNext(const char *opt = "colz");
    void DrawSelected(int, const char *opt = "colz");
    void DrawSearchFrameFunc();
    TString *GetInfo();

    //2---
    void ShowAllButtonFunc(int, const char *opt = "colz");
    void ShowIPButtonFunc(int, const char *opt = "colz");
    void ShowPolButtonFunc(int, const char *opt = "colz");

    //3---
    void NSigPedButtonFunc(double);
    void ShowPedMeanButtonFunc();
    void ShowPedSigmaButtonFunc();
    void AnalysisPed(TString path);

    int GetIp() { return ip; }
    int GetNEvent() { return fEventList.size(); }

    //4--- settings
    void ReadSettings(TGText *text);
    TString GenerateSettingsText();
    TString GenerateSettingsOutput();
    void UpdateEnvParameters(TEnv *env);
    void SaveSettingsToEnv(TEnv *env);

    void SetNEvent(int val) { nEventToAnalysis = val; }
    void DataSwitch(bool flag) { useped = flag; }

private:
    int ip;
    int nEvent;

    //4--- settings parameters
    int nEventToAnalysis;
    int nEtchingMatrix;
    int nExpandMatrix;
    int nEtchExpand;
    int ByMethod;

    bool useped; // flag for using ped or not
    double nped; // nped * sigma will be removed

    //methode-1
    int nXSpatRes;
    int nYSpatRes;
	double nEllipticity;
    double rMinScale;
    double rMaxScale;

    //methode-2
    int HitDist;
    int MinHitsAsCluster;

    TH2F *hPed;
    TH2F *hPedm;
    TH2F *hPeds;

    TString filePath;
    TString rootFilePath;
    TString pedPath;
    TString rootPedPath;

    TH2F *hAll1;       // plot all hits in one histogram
    TH2F *hAll2;       // plot all hits in one histogram (after remove ped)
    TH2F *hBaryCenter; // plot the all barycenter
    TH2F *hIPoint;     // plot the all IP point;
    TH1F *hPol1;       // plot the barycenter line
    TH1F *hPol2;       // plot the polarization

    vector<TH1F *> fPed;
    vector<MyEventClass *> fEventList;
};

#endif
#include "MyRootClass.h"

using namespace std;

MyRootClass::MyRootClass(TString fp, TString pp)
{
    ip = 0;
    hPed = NULL;
    hPedm = NULL;
    hPeds = NULL;
    hAll1 = NULL;
    filePath = fp;
    pedPath = pp;

    nped = 3;
    useped = false;
    nEventToAnalysis = 100;
    nEtchingMatrix = 3;
    nExpandMatrix = 3;
    nEtchExpand = 2;
    ByMethod = 1;
    nXSpatRes = 2;
    nYSpatRes = 2;
    rMinScale = 1.;
    rMaxScale = 1.;
    nEllipticity = 1.5;

    HitDist = 0;
    MinHitsAsCluster = 20;
};

MyRootClass::~MyRootClass()
{
    // Destructor.
}

//______________________________________________________________________________
// analysis settings
//
TString MyRootClass::GenerateSettingsText()
{
    TString settings("#Analysis settings:\n\n");
    settings += TString("#The number of events to be analysised， set to '0' means analysis all events\n");
    settings += TString(Form("nEventToAnalysis = %d\n\n", nEventToAnalysis));
    settings += TString("#Analysis method for reconstruction, should be 1/2, default is 1.\n");
    settings += TString(Form("ByMethod = %d\n\n", ByMethod));
    settings += TString("\n-------------------------------\n");
    settings += TString("--                           --\n");
    settings += TString("--         Methode 1         --\n");
    settings += TString("--   by barycenter line fit  --\n");
    settings += TString("-------------------------------\n\n");
    settings += TString("#The ranking of the etching matrix, should be within [2, 5], default is 3.\n");
    settings += TString(Form("nEtchingMatrix = %d\n\n", nEtchingMatrix));
    settings += TString("#The ranking of the expandng matrix, should be within [2, 5], default is 3.\n");
    settings += TString(Form("nExpandMatrix = %d\n\n", nExpandMatrix));
    settings += TString("#The number of times for etching & expand, shold be positive, default is 2.\n");
    settings += TString(Form("nEtchExpand = %d\n\n", nEtchExpand));
    settings += TString("#The X spatial resolution criteria, 3*sig is used for seperate the 'fatty' events, default is 2 pads.\n");
    settings += TString(Form("nXSpatRes = %d\n\n", nXSpatRes));
    settings += TString("#The Y spatial resolution criteria, 3*sig is used for seperate the 'fatty' events, default is 2 pads.\n");
    settings += TString(Form("nYSpatRes = %d\n\n", nYSpatRes));
    settings += TString("#The Rmin scale value, should be positive > 0, default is 1.\n");
    settings += TString(Form("rMinScale = %1.2f\n\n", rMinScale));
    settings += TString("#The Rmax scale value, should be positive > 0, default is 1.\n");
    settings += TString(Form("rMaxScale = %1.2f\n\n", rMaxScale));
    settings += TString("#The Fatty events screen by ellipticity, should be positive > 1, default is 1.5.\n");
    settings += TString(Form("nEllipticity = %1.2f\n\n", nEllipticity));

    settings += TString("\n-------------------------------\n");
    settings += TString("--                           --\n");
    settings += TString("--         Methode 2         --\n");
    settings += TString("--   by bessel line fit      --\n");
    settings += TString("-------------------------------\n\n");
    settings += TString("#The distance between a cluster and a hit, should be [0~2].\n");
    settings += TString(Form("HitDist = %d\n\n", HitDist));
    settings += TString("#The minimal cluster size, should be lager than 0.\n");
    settings += TString(Form("MinHitsAsCluster = %d\n\n", MinHitsAsCluster));


    return settings;
}

TString WildCardReplace(TString line)
{
    for (int i = 0; i < line.Length(); i++)
    {
        if (line[i] == 46 || (48 <= line[i] && line[i] <= 57))
            continue;
        line[i] = ' ';
    }
    return line;
}

void MyRootClass::ReadSettings(TGText *text)
{
    int nline = 0;
    TGLongPosition pos(0, nline);

    while (text->GetLineLength(nline) != -1)
    {
        if (text->GetChar(pos) != -1)
        {
            TString line(text->GetLine(pos, 100));
            if (line.BeginsWith("nEventToAnalysis"))
                nEventToAnalysis = WildCardReplace(line).Atoi();
            if (line.BeginsWith("nEtchingMatrix"))
                nEtchingMatrix = WildCardReplace(line).Atoi();
            if (line.BeginsWith("nExpandMatrix"))
                nExpandMatrix = WildCardReplace(line).Atoi();
            if (line.BeginsWith("nEtchExpand"))
                nEtchExpand = WildCardReplace(line).Atoi();
            if (line.BeginsWith("ByMethod"))
                ByMethod = WildCardReplace(line).Atoi();
            if (line.BeginsWith("nXSpatRes"))
                nXSpatRes = WildCardReplace(line).Atoi();
            if (line.BeginsWith("nYSpatRes"))
                nYSpatRes = WildCardReplace(line).Atoi();
            if (line.BeginsWith("rMinScale"))
                rMinScale = WildCardReplace(line).Atof();
            if (line.BeginsWith("rMaxScale"))
                rMaxScale = WildCardReplace(line).Atof();
            if (line.BeginsWith("nEllipticity"))
                nEllipticity = WildCardReplace(line).Atof();
            if (line.BeginsWith("HitDist"))
                HitDist = WildCardReplace(line).Atoi();
            if (line.BeginsWith("MinHitsAsCluster"))
                MinHitsAsCluster = WildCardReplace(line).Atoi();
            
        }
        pos.fY = (nline++);
    }
}

TString MyRootClass::GenerateSettingsOutput()
{
    TString settings(Form("Use ped : %s\n", useped ? "YES" : "NO"));
    settings += TString(Form("Pedstal : %.1f * sigma\n", nped));
    settings += TString(Form("nEventToAnalysis : %d\n", nEventToAnalysis));
    settings += TString(Form("nEtchingMatrix : %d\n", nEtchingMatrix));
    settings += TString(Form("nExpandMatrix : %d\n", nExpandMatrix));
    settings += TString(Form("nEtchExpand : %d\n", nEtchExpand));
    settings += TString(Form("ByMethod : %d\n", ByMethod));
    settings += TString(Form("nXSpatRes : %d\n", nXSpatRes));
    settings += TString(Form("nYSpatRes : %d\n", nYSpatRes));
    settings += TString(Form("rMinScale : %1.2f\n", rMinScale));
    settings += TString(Form("rMaxScale : %1.2f\n", rMaxScale));
    settings += TString(Form("nEllipticity : %1.2f\n", nEllipticity));
    settings += TString(Form("HitDist : %d\n", HitDist));
    settings += TString(Form("MinHitsAsCluster: %d\n", MinHitsAsCluster));

    return settings;
}

void MyRootClass::UpdateEnvParameters(TEnv *env)
{
    nEventToAnalysis = env->GetValue("nEventToAnalysis", nEventToAnalysis);
    nEtchingMatrix = env->GetValue("nEtchingMatrix", nEtchingMatrix);
    nExpandMatrix = env->GetValue("nExpandMatrix", nExpandMatrix);
    nEtchExpand = env->GetValue("nEtchExpand", nEtchExpand);
    ByMethod = env->GetValue("ByMethod", ByMethod);
    nXSpatRes = env->GetValue("nXSpatRes", nXSpatRes);
    nYSpatRes = env->GetValue("nYSpatRes", nYSpatRes);
    rMinScale = env->GetValue("rMinScale", rMinScale);
    rMaxScale = env->GetValue("rMaxScale", rMaxScale);
    nEllipticity = env->GetValue("nEllipticity", nEllipticity);
    HitDist = env->GetValue("HitDist", HitDist);
    MinHitsAsCluster = env->GetValue("MinHitsAsCluster", MinHitsAsCluster);
}

void MyRootClass::SaveSettingsToEnv(TEnv *env)
{
    env->SetValue("nEventToAnalysis", nEventToAnalysis);
    env->SetValue("nEtchingMatrix", nEtchingMatrix);
    env->SetValue("nExpandMatrix", nExpandMatrix);
    env->SetValue("nEtchExpand", nEtchExpand);
    env->SetValue("ByMethod", ByMethod);
    env->SetValue("nXSpatRes", nXSpatRes);
    env->SetValue("nYSpatRes", nYSpatRes);
    env->SetValue("rMinScale", rMinScale);
    env->SetValue("rMaxScale", rMaxScale);
    env->SetValue("nEllipticity", nEllipticity);
    env->SetValue("HitDist", HitDist);
    env->SetValue("MinHitsAsCluster", MinHitsAsCluster);
    env->SaveLevel(kEnvLocal);
}

//______________________________________________________________________________
// analysis pedstal
//
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
            for (int k = 0; k < (int)dvec[i][j].size(); k++)
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
            for (int k = 0; k < (int)dvec[i][j].size(); k++)
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
        for (int j = 0; j < NY; j++)
            fPed[i * NY + j]->Write();

    cout << "---->Save ped root file to: " << rootPedPath << endl;
    rootPedFile->Close();
    hPedm->Draw("colz");
}

//______________________________________________________________________________
// use nsig * pedstal_sigma as threshold
//
void MyRootClass::NSigPedButtonFunc(double _ped)
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

void MyRootClass::ShowPedMeanButtonFunc()
{
    if (hPedm != NULL)
        hPedm->Draw("colz");
}

void MyRootClass::ShowPedSigmaButtonFunc()
{
    if (hPeds != NULL)
        hPeds->Draw("colz");
}

//______________________________________________________________________________
// read ped distribution from data/root file
//
void MyRootClass::ReadPed()
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
        if (rootPedFile->IsOpen())
        {
            hPedm = (TH2F *)rootPedFile->Get("hPedm");
            hPeds = (TH2F *)rootPedFile->Get("hPeds");
            hPed = (TH2F *)rootPedFile->Get("hPeds");

            for (int i = 0; i < NX; i++)
                for (int j = 0; j < NY; j++)
                    hPed->SetBinContent(i + 1, j + 1, hPedm->GetBinContent(i + 1, j + 1) + nped * hPeds->GetBinContent(i + 1, j + 1));
        }
        else
        {
            hPed = NULL;
            hPedm = NULL;
            hPeds = NULL;
        }
    }
}

//______________________________________________________________________________
// init histogram for plot
//
void MyRootClass::IniHistogram()
{
    ip = 0;
    nEvent = 0;

    if (hAll1 != NULL)
    {
        delete hAll1;
        delete hAll2;
        delete hBaryCenter;
        delete hIPoint;
        delete hPol1;
        delete hPol2;
    }
    hAll1 = new TH2F("hAll1", "Histogram of all hits", NX, 1, NX + 1, NY, 1, NY + 1);
    hAll2 = new TH2F("hAll2", "Histogram of all hits", NX, 1, NX + 1, NY, 1, NY + 1);
    hBaryCenter = new TH2F("hBaryCenter", "Histogram of all barycenters", NX, 1, NX + 1, NY, 1, NY + 1);
    hIPoint = new TH2F("hIPoint", "Histogram of all IP points", NX, 1, NX + 1, NY, 1, NY + 1);
    hPol1 = new TH1F("hPol1", "#theta of the bary-center line", 100, -TMath::Pi(), TMath::Pi());
    hPol2 = new TH1F("hPol2", "#theta of the fitted Polarization \n par[0]*pow(cos(x[0]-par[1]),2)+par[2]", 40, -TMath::Pi(), TMath::Pi());
    hPol2->GetXaxis()->SetTitle("#theta");
    hPol2->GetYaxis()->SetTitle("counts");
    fEventList.clear();
}

//______________________________________________________________________________
// analysis directory
//
void MyRootClass::AnalysisDir(TString fileDir)
{
    if (useped)
        ReadPed();
    IniHistogram();

    FILE *fp = gSystem->OpenPipe("ls " + fileDir + "/*.data", "r");
    if (!fp)
    {
        cout << "----> NO data files exists in " << fileDir << "!" << endl;
        return;
    }

    vector<TString> dataList;
    char line[1000];
    while (fgets(line, sizeof(line), fp))
    {
        TString s(line);
        if (s.Index(".data") == -1)
            continue;
        dataList.push_back(s.ReplaceAll("\n", ""));
    }
    cout << "----> " << dataList.size() << " data files exist in " << fileDir << "." << endl;

    for (int i = 0; i < (int)dataList.size(); i++)
        AnalysisFile(dataList[i]);

    cout << "----> " << nEvent << " events recorded." << endl;
}

//______________________________________________________________________________
// analysis a single file
//
void MyRootClass::Analysis(TString filePath)
{
    if (useped)
        ReadPed();
    IniHistogram();
    AnalysisFile(filePath);

    cout << "----> " << nEvent << " events recorded." << endl;
}

//______________________________________________________________________________
// analysis a single file
//
void MyRootClass::AnalysisFile(TString filePath)
{
    cout << "--> Opening: " << filePath << endl;
    ifstream ifSignal(filePath, ios::binary);

    while (ifSignal.good())
    {
        unsigned short _data[NX][NY];
        ifSignal.read((char *)(&_data), sizeof(_data));

        if (_data[0][0] == 0xFFFF && _data[1][1] == 0xFFFF)
        {
            ifSignal.read((char *)(&_data), sizeof(_data));
            continue;
        }

        MyEventClass *fEvent = new MyEventClass(nEvent, 0, NX - 1, 0, NY - 1);
        
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
            {
                fEvent->Fill(i, j, _data[i][j]);
                double q = _data[i][j]; // - hPed->GetBinContent(i + 1, j + 1);
                hAll1->Fill(i + 1, j + 1, (q < 0) ? 0 : q);
            }

        fEvent->SetEtchingMatrixRank(nEtchingMatrix);
        fEvent->SetExpandMatrixRank(nExpandMatrix);
        fEvent->SetNEtchExpand(nEtchExpand);
        fEvent->SetMethod(ByMethod);

        //设置分析用的参数控制
        if (ByMethod == 2)
        {
            fEvent->SetHitDist(HitDist);
            fEvent->SetMinHitsAsCluster(MinHitsAsCluster);
        }
        else //method 1
        {
            fEvent->SetFattyEventCutOnX(3 * nXSpatRes);
            fEvent->SetFattyEventCutOnY(3 * nYSpatRes);
            fEvent->SetRMinScale(rMinScale);
            fEvent->SetRMaxScale(rMaxScale);
            fEvent->SetEllipticity(nEllipticity);
        }

        //分析
        fEvent->GenerateHist((useped) ? hPed : NULL);

        //分析结果填图
        if (ByMethod == 2)
        {
			fEvent->FillPolarization2(hPol2);
        }
        else //method 1
        {
            fEvent->Fill2DPlot(hAll2);
            fEvent->FillBaryCenter(hBaryCenter);
            fEvent->FillIPpoint(hIPoint);
			fEvent->FillPolarization1(hPol1);
            fEvent->FillPolarization2(hPol2);
        }

        fEventList.push_back(fEvent);

        nEvent++;
        //delete fEvent;

        if (nEvent % 100 == 0)
            cout << "Dealed: " << nEvent << " " << nEventToAnalysis << endl;

        if (nEventToAnalysis > 0 && nEvent >= nEventToAnalysis)
            break;
    }

    ifSignal.close();
}

//______________________________________________________________________________
// draw events
//
void MyRootClass::DrawRaw(const char *opt)
{
    fEventList.at(ip)->Get2DRawPlot()->Draw(opt);
}

void MyRootClass::DrawPre(const char *opt)
{
    if (ip <= 0)
        return;

    int _ip = ip;
    while (_ip >= 1)
    {
        _ip--;
        if (fEventList.at(_ip)->GetDataQuality() == 1)
            break;
    }
    ip = _ip;
    fEventList.at(ip)->Draw2DResult(opt);
}

void MyRootClass::DrawNext(const char *opt)
{
    if ((unsigned long)ip == fEventList.size())
        return;

    int _ip = ip;
    while (_ip < (int)fEventList.size() - 1)
    {
        _ip++;
        //if (fEventList.at(_ip)->GetNeedNum() < 0.1 && fEventList.at(_ip)->GetNeedNum() > 0){  //check the special dot
        //cout<<fEventList.at(_ip)->GetNeedNum()<<endl;
        if (fEventList.at(_ip)->GetDataQuality() == 1)
            break;
        //}
    }

    ip = _ip;
    fEventList.at(ip)->Draw2DResult(opt);
}

void MyRootClass::DrawSelected(int _ip, const char *opt)
{
    if (_ip < 0 || _ip > (int)fEventList.size() - 1)
        return;

    ip = _ip;
    fEventList.at(ip)->Draw2DResult(opt);
}

TString *MyRootClass::GetInfo()
{
    return fEventList.at(ip)->GetInfo();
}

void MyRootClass::DrawSearchFrameFunc()
{
    fEventList.at(ip)->DrawSearchResult();
}

//______________________________________________________________________________
// show events
//
void MyRootClass::ShowAllButtonFunc(int flag, const char *opt)
{
    if (flag == 1)
        hAll1->Draw(opt);
    else
        hAll2->Draw(opt);
}

void MyRootClass::ShowIPButtonFunc(int flag, const char *opt)
{
    if (flag == 1)
        hBaryCenter->Draw(opt);
    else
        hIPoint->Draw(opt);
}

//______________________________________________________________________________
// fit the polarization distribution
//
Double_t fitFunction(Double_t *x, Double_t *par)
{
    //	return abs(par[0]*sin(par[1]*x[0]))+par[2];
    return par[0] * pow(cos(x[0] - par[1]), 2) + par[2];
}

void MyRootClass::ShowPolButtonFunc(int flag, const char *opt)
{
    TF1 *hisFunc = new TF1("hisFunc", fitFunction, -TMath::Pi(), TMath::Pi(), 3);
    //TF1 *hisFunc = new TF1("hisFunc",fitLineFunction,0,TMath::Pi(),2);
    ofstream outf;
    outf.open("./abc.txt", ios::app);
    if (flag == 1)
        hPol1->Draw(opt);
    else
    {
        hisFunc->SetParameters(100, 1, 100);
        //hisFunc->SetParameters(0,200);
        hPol2->Fit(hisFunc, "r");
        hPol2->SetMinimum(200);
        hPol2->Draw(opt);
        hPol2->SetStats(kTRUE);
        gStyle->SetOptFit(1);
        double param1 = hisFunc->GetParameter(0); //输出到txt文件
        outf << param1 << " ";
        double param2 = hisFunc->GetParameter(2);
        outf << param2 << " ";
        double parError1 = hisFunc->GetParError(0);
        outf << parError1 << " ";
        double parError2 = hisFunc->GetParError(2);
        outf << parError2 << endl;
        //double thechi2=hisFunc->GetChisquare();
        //outf<<thechi2<<" "<<endl;
        //cout<<"_____ "<<param1/(2*param2+param1)<<endl;
    }
}

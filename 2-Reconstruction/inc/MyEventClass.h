#ifndef _MyEventClass_h_
#define _MyEventClass_h_

#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TList.h"
#include "TLine.h"
#include "TMarker.h"
#include "TGaxis.h"
#include "TEllipse.h"
#include "TFile.h"

#include <vector>
#include <map>
#include <iostream>

//Quality of data, set according to criteria
#define QF_PED 0
#define QF_HIT 1
#define QF_FAT 2
#define QF_MuLike 3
#define QF_AfterPulse 4

using namespace std;

//-----------------------------------------------------------------------------------------------
//
class MyEventClass
{
public:
    MyEventClass(int _id, int _xmin, int _xmax, int _ymin, int _ymax);
    virtual ~MyEventClass();

    void Fill(int x, int y, double _b) { data[x - xmin][y - ymin] = _b; };
    void GenerateHist(TH2F *ped, bool anaflag = kTRUE);

    vector<vector<int>> GetData() { return data; }
    double GetData(int x, int y) { return data[x - xmin][y - ymin]; }

    //--------------------------
    // 算法1:
    void AnalysisHist();
    void Fill2DPlot(TH2F *);
    void FillBaryCenter(TH2F *);
    void FillIPpoint(TH2F *);
    void FillPolarization1(TH1F *);
    void FillPolarization2(TH1F *);
    void Filllengththelongest(TH1F *);

    TH2F *Get2DPlot() { return f2D; };
    TH2F *Get2DRawPlot() { return f2D_raw; };
    TMarker *GetBaryCenterAsMarker() { return mBcenter; }
    TLine *GetPrincipalAxis1() { return lPrinAxis1; }
    TLine *GetPrincipalAxis2() { return lPrinAxis2; }
    TMarker *GetConvertionPoint() { return mCovPoint; }
    TLine *GetCovertionAxis() { return lCovAxis; }
    TString *GetInfo() { return info; };

    void SetDataQuality(int f) { dataQFlag = f; }
    int GetDataQuality() { return dataQFlag; }
    int GetClusterSize() { return clusterSize; }
    int GetPulseHeight() { return pulseHeight; }
    double GetNeedNum() { return aTheta1; } // can choose a value

    //--- 算法1的相关参数设置
    void SetEtchingMatrixRank(int rank) { nEtchingMatrix = rank; }
    void SetExpandMatrixRank(int rank) { nExpandMatrix = rank; }
    void SetNEtchExpand(int times) { nEtchExpand = times; }
    void SetFattyEventCutOnX(int xf) { nFattyX = xf; }
    void SetFattyEventCutOnY(int yf) { nFattyY = yf; }
    void SetEllipticity(double nEllip) { nEllipticity = nEllip; }
    void SetRMinScale(double scale) { rMinScale = scale; }
    void SetRMaxScale(double scale) { rMaxScale = scale; }

    //--- 算法1的作图
    void Draw2DResultMethod1(const char *opt);
    void DrawSearchResultMethod1();

    //--------------------------
    //算法2
    void AnalysisHist2();

    void Draw2DResultMethod2(const char *opt);
    void DrawSearchResultMethod2() {;}

    //--------------------------
    //结果
    void Draw2DResult(const char *opt) { Draw2DResultMethod2(opt); }
    void DrawSearchResult() { DrawSearchResultMethod2(); }

private:
    //--- basic info.
    int id;
    int nx, ny;
    int xmin, xmax;
    int ymin, ymax;
    int clusterSize;
    double pulseHeight;

    int dataQFlag; // data quality
    vector<vector<int>> data;

    //--- analysis var.
    //算法1
    //double mBx, mBy;
    double mCx, mCy;
    double lCk, lCb;
    double aTheta1, aTheta2;
    double mom2nd;
    double mom2ndMax;
    double mom2ndMin;
    double defVal;

    TMarker *mBcenter;
    TMarker *mCovPoint;
    TLine *lPrinAxis1;
    TLine *lPrinAxis2;
    TEllipse *e1;
    TEllipse *e2;
    TLine *lCovAxis;

    //算法2
    vector<pair<double, double>> plist;

    //--- 通用.
    TH2F *f2D;
    TH2F *f2D_raw;

    TString *info;

    //---- etch&expand
    int nEtchingMatrix;
    int nExpandMatrix;
    int nEtchExpand;
    int nFattyX;
    int nFattyY;
    double nEllipticity;
    double rMinScale;
    double rMaxScale;

    void EtchHistogram(TH2F *, TH2F *);
    void ExpandHistogram(TH2F *, TH2F *);
};

#endif
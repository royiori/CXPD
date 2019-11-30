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
#include "Riostream.h"

#include <vector>
#include <map>

//Quality of data, set according to criteria
#define QF_PED 0
#define QF_HIT 1
#define QF_FAT 2
#define QF_MultiCluster 3
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

    void SetMethod(int m) { ByMethod = m; }
    TH2F *Get2DPlot() { return f2D; };
    TH2F *Get2DRawPlot() { return f2D_raw; };

    //--------------------------
    // 算法1:
    void AnalysisHist1();
    void Fill2DPlot(TH2F *);
	void FillBaryCenter(TH2F *h) { if (mBx != defVal && mBy != defVal) h->Fill(mBx, mBy); }
    void FillIPpoint(TH2F *h) {if (mBx2 != defVal && mBy2 != defVal) h->Fill(mBx2, mBy2); }
    void FillPolarization1(TH1F *h) { if (aTheta1 != defVal) h->Fill(aTheta1); }
    void FillPolarization2(TH1F *h) { if (aTheta2 != defVal) {
		//if (aTheta2 < 0)
		//aTheta2 = aTheta2 + TMath::Pi();
		h->Fill(aTheta2); }}
    void Filllengththelongest(TH1F *h) { h->Fill(clusterlength); }

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
    double GetNeedNum() { return llm; } // can choose a value

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
    vector<pair<double, double>> AnalysisCluster2(vector<pair<double, double>> hlist);

    //--- 算法2的相关参数设置
    void SetHitDist(int hd) { HitDist = hd; }
    void SetMinHitsAsCluster(int nh) { MinHitsAsCluster = nh; }

    //--- 算法2的作图
    void Draw2DResultMethod2(const char *opt);
    void DrawSearchResultMethod2();

    //--------------------------
    //结果
    void Draw2DResult(const char *opt)
    {
        if (ByMethod == 2)
            Draw2DResultMethod2(opt);
        else
            Draw2DResultMethod1(opt);
    }
    void DrawSearchResult()
    {
        if (ByMethod == 2)
            DrawSearchResultMethod2();
        else
            DrawSearchResultMethod1();
    }

private:
    //--- basic info.
    int id;
    int nx, ny;
	int llm;
    int xmin, xmax;
    int ymin, ymax;
    int clusterSize;
    double pulseHeight;

    int dataQFlag; // data quality
    vector<vector<int>> data;

    int ByMethod;
    
    //--- analysis var.
    //算法1
	double mBx, mBy, mBx2, mBy2;
    //double mCx, mCy;
    double lCk, lCb;
    double aTheta1, aTheta2;
    double mom2nd;
    double mom2ndMax;
    double mom2ndMin;
    double defVal;
	double clusterlength;
	double pureMom2nd;
	double xn;

    TMarker *mBcenter;
    TMarker *mCovPoint;
    TLine *lPrinAxis1;
    TLine *lPrinAxis2;
	TLine *lP0p1;
	TLine *lP1p2;
	TLine *lP2p3;
	TLine *lPrinAxisN[30];
    TEllipse *e1;
    TEllipse *e2;
    TLine *lCovAxis;

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

    //算法2
    vector<vector<pair<double, double>>> clist; //数据分cluster后的列表
    vector<vector<pair<double, double>>> plist; //Bessel拟合后的plist数据
    int HitDist;
    int MinHitsAsCluster;

    //--- 通用.
    TH2F *f2D_raw; //原始数据
    TH2F *f2D; //算法1是etch以后的数据，算法2是找到cluster以后的数据

    TString *info;
};

#endif

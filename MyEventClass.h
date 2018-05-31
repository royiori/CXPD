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

#include <vector>
#include <map>
#include <iostream>

//----
// fit data
//
double xmean = 0; //bary
double ymean = 0;
double qtot = 0;
double xmean2 = 0; //impact point
double ymean2 = 0;
double qtot2 = 0;

const int QF_Ped = 0;
const int QF_Hit = 1;
const int QF_MuLike = 2;
const int QF_AfterPulse = 3;

//----
// chisquare function
//
std::vector<std::pair<double, double>> coords;
std::vector<double> values;
std::vector<double> errors;

void myFcn(Int_t & /*nPar*/, Double_t * /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)
{
    int n = coords.size();
    double chi2 = 0;
    double tmp, x, y;
    double dist;

    double k = p[0];
    double b = ymean - k * xmean;

    for (int i = 0; i < n; ++i)
    {
        if (values[i] == 0)
            continue;
        if (p[0] == 0 && p[1] == 0)
            continue;

        x = coords[i].first;
        y = coords[i].second;

        dist = fabs(k * x + b - y) / sqrt(1 + k * k);

        chi2 += pow(values[i] / qtot * dist, 2);
    }
    fval = (chi2 == 0) ? 1E19 : chi2;
}

void myFcn2(Int_t & /*nPar*/, Double_t * /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)
{
    int n = coords.size();
    double chi2 = 0;
    double tmp, x, y;
    double dist;

    double k = p[0];
    double b = ymean2 - k * xmean2;

    for (int i = 0; i < n; ++i)
    {
        if (values[i] == 0)
            continue;
        if (p[0] == 0 && p[1] == 0)
            continue;

        x = coords[i].first;
        y = coords[i].second;

        dist = fabs(k * x + b - y) / sqrt(1 + k * k);

        double d0 = sqrt((x - xmean2) * (x - xmean2) + (y - ymean2) * (y - ymean2));
        //d0 = (d0<1) ? 1 : d0;

        chi2 += pow(dist * exp(-d0/2), 2);
    }
    fval = (chi2 == 0) ? 1E19 : chi2;
}

//-----------------------------------------------------------------------------------------------
//
class MyEventClass
{
  public:
    MyEventClass(int _id, int _xmin, int _xmax, int _ymin, int _ymax)
    {
        dataQFlag = QF_Ped;
        id = _id;
        xmin = _xmin;
        xmax = _xmax;
        ymin = _ymin;
        ymax = _ymax;
        nx = xmax - xmin + 1;
        ny = ymax - ymin + 1;
        f2D = NULL;
        f2D_raw = NULL;
        mBcenter = NULL;
        mCovPoint = NULL;
        lPrinAxis1 = NULL;
        lPrinAxis2 = NULL;
        e1 = NULL;
        e2 = NULL;
        lCovAxis = NULL;
        info = NULL;

        data.resize((nx));
        for (int i = 0; i < nx; i++)
            data[i].resize(ny);
    }

    virtual ~MyEventClass();

    void Fill(int x, int y, double _b) { data[x - xmin][y - ymin] = _b; };
    void GenerateHist(TH2F *ped, bool anaflag);
    void AnalysisHist();
    void Fill2DPlot(TH2F *);
    void FillBaryCenter(TH2F *);
    void FillIPpoint(TH2F *);


    TH2F *Get2DPlot() { return f2D; };
    TH2F *Get2DRawPlot() { return f2D_raw; };

    TMarker *GetBaryCenterAsMarker() { return mBcenter; }
    TLine *GetPrincipalAxis1() { return lPrinAxis1; }
    TLine *GetPrincipalAxis2() { return lPrinAxis2; }
    TMarker *GetConvertionPoint() { return mCovPoint; }
    TLine *GetCovertionAxis() { return lCovAxis; }
    TString *GetInfo() { return info; };

    void SetDataQuality(int f) {dataQFlag = f;}
    int GetDataQuality() { return dataQFlag; }
    int GetClusterSize() { return clusterSize;}
    int GetPulseHeight() { return pulseHeight;}
    vector<vector<int>> GetData() {return data;}
    double GetData(int x, int y) { return data[x-xmin][y-ymin]; }

    void Draw2DResult(const char *opt)
    {
        if(f2D!=NULL) f2D->Draw(opt);
        if(mBcenter!=NULL) mBcenter->Draw();
        if(mCovPoint!=NULL) mCovPoint->Draw();
        if(lPrinAxis1!=NULL) lPrinAxis1->Draw();
        if(lCovAxis!=NULL) lCovAxis->Draw();
    }

    void DrawEllipse()
    {
        if(lPrinAxis2!=NULL) lPrinAxis2->Draw();
        if(e1!=NULL) e1->Draw();
        if(e2!=NULL) e2->Draw();
    }

  private:
    int xmin;
    int xmax;
    int ymin;
    int ymax;
    int nx;
    int ny;
    int id;

    int clusterSize;
    int dataQFlag; // QF_Ped/QF_Hit/...
    double pulseHeight;

    double mBx, mBy;
    double lPk, lPb;
    double mCx, mCy;
    double lCk, lCb;
    double mom2nd;
    double mom3rd;

    TH2F *f2D;
    TH2F *f2D_raw;

    TMarker *mBcenter;
    TMarker *mCovPoint;
    TLine *lPrinAxis1;
    TLine *lPrinAxis2;
    TEllipse *e1;
    TEllipse *e2;
    TLine *lCovAxis;
    TString *info;

    vector<vector<int>> data;

    void EtchHistogram(TH2F *, TH2F *);
    void ExpandHistogram(TH2F *, TH2F *);
};

MyEventClass::~MyEventClass()
{
    // Destructor.
    if(f2D       !=NULL) delete f2D        ;
    if(f2D_raw   !=NULL) delete f2D_raw    ;
    if(mBcenter  !=NULL) delete mBcenter   ;
    if(mCovPoint !=NULL) delete mCovPoint  ;
    if(lPrinAxis1!=NULL) delete lPrinAxis1 ;
    if(lPrinAxis2!=NULL) delete lPrinAxis2 ;
    if(e1        !=NULL) delete e1         ;
    if(e2        !=NULL) delete e2         ;
    if(lCovAxis  !=NULL) delete lCovAxis   ;
    if(info      !=NULL) delete info       ;

    data.clear();
}

//______________________________________________________________________________
void MyEventClass::EtchHistogram(TH2F *f0, TH2F *f1)
{
    const int netch = 3;
    int ietch = 1;
    int etchData[netch][netch] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    //... etching
    for (int i = xmin; i < nx; i++)
        for (int j = ymin; j < ny; j++)
            f1->SetBinContent(i + 1, j + 1, 0);

    for (int i = xmin; i < nx; i++)
        for (int j = ymin; j < ny; j++)
        {
            int imin = (i - ietch < 0) ? 0 : i - ietch;
            int jmin = (j - ietch < 0) ? 0 : j - ietch;
            int imax = (i + ietch > nx - ietch) ? nx - ietch : i + ietch;
            int jmax = (j + ietch > ny - ietch) ? ny - ietch : j + ietch;

            int flag = 1;
            for (int ii = imin; ii <= imax; ii++)
                for (int jj = jmin; jj <= jmax; jj++)
                {
                    int m = ii - imin;
                    int n = jj - jmin;
                    if (etchData[m][n] == 0)
                        continue;
                    flag *= etchData[m][n] * f0->GetBinContent(ii + 1, jj + 1);
                }

            if (flag == 0)
                f1->SetBinContent(i + 1, j + 1, 0);
            else
                f1->SetBinContent(i + 1, j + 1, 1);
        }
    return;
}

void MyEventClass::ExpandHistogram(TH2F *f0, TH2F *f1)
{
    const int netch = 3;
    int ietch = 1;
    int etchData[netch][netch] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    //... expand
    for (int i = xmin; i < nx; i++)
        for (int j = ymin; j < ny; j++)
            f1->SetBinContent(i + 1, j + 1, 0);

    for (int i = xmin; i < nx; i++)
        for (int j = ymin; j < ny; j++)
        {
            if (f0->GetBinContent(i + 1, j + 1) == 0)
                continue;

            int imin = (i - ietch < 0) ? 0 : i - ietch;
            int jmin = (j - ietch < 0) ? 0 : j - ietch;
            int imax = (i + ietch > nx - ietch) ? nx - ietch : i + ietch;
            int jmax = (j + ietch > ny - ietch) ? ny - ietch : j + ietch;

            for (int ii = imin; ii <= imax; ii++)
                for (int jj = jmin; jj <= jmax; jj++)
                {
                    int m = ii - imin;
                    int n = jj - jmin;
                    if (etchData[m][n] == 0)
                        continue;
                    f1->SetBinContent(ii + 1, jj + 1, 1);
                }
        }
    return;
}

//______________________________________________________________________________
void MyEventClass::GenerateHist(TH2F *hPed, bool anaflag = kTRUE)
{
    if (f2D != NULL)
    {
        delete f2D;
        delete f2D_raw;
    }

    if (gDirectory->Get(Form("f2D_%d", id)) != NULL)
    {
        TH2F *f1 = (TH2F *)gDirectory->Get(Form("f2D_%d", id));
        TH2F *f2 = (TH2F *)gDirectory->Get(Form("f2D_%d_0", id));
        delete f1;
        delete f2;
    }

    f2D = new TH2F(Form("f2D_%d", id), Form("Event %d", id), nx, 1, xmax - xmin + 1, ny, 1, ymax - ymin + 1);
    f2D_raw = new TH2F(Form("f2D_%d_0", id), Form("Event %d", id), nx, 1, xmax - xmin + 1, ny, 1, ymax - ymin + 1);

    //... store to data
    coords.clear();
    values.clear();
    errors.clear();
    for (int i = xmin; i < nx; i++)
        for (int j = ymin; j < ny; j++)
        {
            double ped = 0;
            if (hPed != NULL)
                ped = hPed->GetBinContent(i + 1, j + 1);

            double q = data[i - xmin][j - ymin] - ped;
            q = (q < 0) ? 0 : q;
            f2D_raw->Fill(i + 1, j + 1, q);
        }

    //... deal histogram
    TH2F *ftmp0 = new TH2F("ftmp0", "ftmp", nx, 1, xmax - xmin + 1, ny, 1, ymax - ymin + 1);
    TH2F *ftmp1 = new TH2F("ftmp1", "ftmp", nx, 1, xmax - xmin + 1, ny, 1, ymax - ymin + 1);
    EtchHistogram(f2D_raw, ftmp0);
    ExpandHistogram(ftmp0, ftmp1);
    ExpandHistogram(ftmp1, ftmp0);
    EtchHistogram(ftmp0, ftmp1);

    //... store to data
    clusterSize = 0;
    pulseHeight = 0;

    for (int i = xmin; i < nx; i++)
        for (int j = ymin; j < ny; j++)
            if (ftmp1->GetBinContent(i + 1, j + 1) > 0)
            {
                clusterSize++;
                pulseHeight += f2D_raw->GetBinContent(i + 1, j + 1);
                f2D->SetBinContent(i + 1, j + 1, f2D_raw->GetBinContent(i + 1, j + 1));

                //.. stored for fit..
                coords.push_back(std::make_pair(i + 1, j + 1));
                values.push_back( f2D_raw->GetBinContent(i + 1, j + 1) );
                errors.push_back(sqrt( f2D_raw->GetBinContent(i + 1, j + 1) ));
            }

    delete ftmp0;
    delete ftmp1;

    //... generate info
    info = new TString();
    info->Append(Form("Event Number:      \t%d\n", id));
    info->Append(Form("Cluster Size:         \t%d\n", clusterSize));
    info->Append(Form("Pulse Height:      \t%.2f\n", pulseHeight));

    dataQFlag = QF_Ped;
    if (hPed==NULL || clusterSize < 20 || pulseHeight < 200)
        return;
    dataQFlag = QF_Hit;
    if(anaflag) AnalysisHist();
    return;
}

//______________________________________________________________________________
void MyEventClass::AnalysisHist()
{
    //----------------------
    //1. find barycenter
    xmean = 0;
    ymean = 0;
    qtot = 0;

    int N = coords.size();
    for (int i = 0; i < N; i++)
    {
        double x = coords[i].first;
        double y = coords[i].second;
        double q = values[i];

        //cout<<i<<": "<<x<<", "<<y<<"  - "<<q<<endl;

        qtot += q;
        xmean += q * x;
        ymean += q * y;
    }

    xmean /= qtot;
    ymean /= qtot;
    mBx = xmean;
    mBy = ymean;

    //1.1 draw result
    mBcenter = new TMarker(mBx, mBy, 30);
    mBcenter->SetMarkerColor(kRed);
    mBcenter->SetMarkerSize(2.6);
    mBcenter->Draw();

    //----------------------
    //2. fit barycenter-line
    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter *minuit = TVirtualFitter::Fitter(0, 10);
    minuit->SetParameter(0, "k", 0, 0.01, 0, 0);
    minuit->SetFCN(myFcn);

    double arglist[100];
    arglist[0] = 0;
    //minuit->ExecuteCommand("SET PRINT", arglist, 2);

    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.01; // tolerance
    minuit->ExecuteCommand("MIGRAD", arglist, 2);

    //2.1 get result
    double chi2, edm, errdef;
    int nvpar, nparx;
    minuit->GetStats(chi2, edm, errdef, nvpar, nparx);

    double k = minuit->GetParameter(0);
    double b = ymean - k * xmean;
    lPk = k;
    lPk = b;

    //2.2 draw result
    double lx1 = (k == 0) ? xmin : (ymin - b) / k;
    double ly1 = (k == 0) ? ymean : ymin;
    double lx2 = (k == 0) ? xmax : (ymax - b) / k;
    double ly2 = (k == 0) ? ymean : ymax;
    double theta1 = atan(k);

    lPrinAxis1 = new TLine(lx1, ly1, lx2, ly2);
    lPrinAxis1->SetLineColor(kBlack);
    lPrinAxis1->Draw();

    //----------------------
    //3. cal mom3rd
    double k2 = (k == 0) ? 1 : -1 / k;
    double b2 = (k == 0) ? 1 : ymean - k2 * xmean;
    double theta2 = (k == 0) ? TMath::Pi()/2 : atan(-1/k);

    mom3rd = 0;
    for (int i = 0; i < N; i++)
    {
        double x = coords[i].first;
        double y = coords[i].second;
        double q = values[i];
        double d = (k == 0) ? (x - xmean) * q / qtot : (y - b2 - k2 * x) / sqrt(1 + k2 * k2) * q / qtot;
        double d3 = d * d * d;
        mom3rd += d3;
    }

    //3.1 draw result
    lx1 = (k == 0) ? xmean : (ymin - b2) / k2;
    ly1 = ymin;
    lx2 = (k == 0) ? xmean : (ymax - b2) / k2;
    ly2 = ymax;
    lPrinAxis2 = new TLine(lx1, ly1, lx2, ly2);
    lPrinAxis2->SetLineColor(kBlue);
    lPrinAxis2->Draw();

    //3.2 cal the average distance
    double dmax = 0;
    int nd = 0;
    for (int i = 0; i < N; i++)
    {
        double x = coords[i].first;
        double y = coords[i].second;
        double q = values[i];
        double d = (k == 0) ? (x - xmean) : (y - b2 - k2 * x) / sqrt(1 + k2 * k2);
        double d3 = d * d * d;
        if (d3 * mom3rd > 0)
            continue;
        dmax += d;
        nd++;
    }
    dmax /= nd;
    dmax *= 1.1;

    //3.3 draw TEllipse
    double phimin = theta2 / TMath::Pi() * 180;
    double phimax = theta2 / TMath::Pi() * 180;
    if(k2>0 && mom3rd>0) phimin -= 180;
    if(k2>0 && mom3rd<0) phimax += 180;
    if(k2<0 && mom3rd>0) phimin -= 180;
    if(k2<0 && mom3rd<0) phimax += 180;

    e1 = new TEllipse(xmean, ymean, dmax, dmax, phimin, phimax);
    e2 = new TEllipse(xmean, ymean, dmax * 2.33, dmax * 2.33, phimin, phimax);
    e1->SetFillColorAlpha(kWhite, 0);
    e2->SetFillColorAlpha(kWhite, 0);
    e1->Draw();
    e2->Draw();

    //----------------------
    //4. find IP
    xmean2 = 0;
    ymean2 = 0;
    mom2nd = 0;
    qtot2  = 0;

    for (int i = 0; i < N; i++)
    {
        double x = coords[i].first;
        double y = coords[i].second;
        double q = values[i];
        double d = (k == 0) ? (x - xmean) : (y - b2 - k2 * x) / sqrt(1 + k2 * k2);
        double d3 = d * d * d;
        if (d3 * mom3rd > 0)
            continue;
        if (sqrt((x - xmean) * (x - xmean) + (y - ymean) * (y - ymean)) < dmax)
            continue;
        xmean2 += q * x;
        ymean2 += q * y;
        qtot2 += q;
    }
    xmean2 /= qtot2;
    ymean2 /= qtot2;
    mCx = xmean2;
    mCy = ymean2;

    mCovPoint = new TMarker(mCx, mCy, 31);
    mCovPoint->SetMarkerColor(kRed);
    mCovPoint->SetMarkerSize(2.6);
    mCovPoint->Draw();

    //----------------------
    //5. Fit eject line
    TVirtualFitter *minuit2 = TVirtualFitter::Fitter(0, 10);
    minuit2->SetParameter(0, "k", 0, 0.01, 0, 0);
    minuit2->SetFCN(myFcn2);

    arglist[0] = 0;
    minuit2->ExecuteCommand("SET PRINT", arglist, 2);

    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.01; // tolerance
    minuit2->ExecuteCommand("MIGRAD", arglist, 2);

    //5.1 get result
    minuit2->GetStats(chi2, edm, errdef, nvpar, nparx);

    double k3 = minuit2->GetParameter(0);
    double b3 = ymean2 - k3 * xmean2;

    lCk = k3;
    lCb = b3;

    //5.2 draw result
    lx1 = (k3 == 0) ? xmean2 : (ymin - b3) / k3;
    ly1 = ymin;
    lx2 = (k3 == 0) ? xmean2 : (ymax - b3) / k3;
    ly2 = ymax;

    lCovAxis = new TLine(lx1, ly1, lx2, ly2);
    lCovAxis->SetLineColor(kRed);
    lCovAxis->Draw();

    //-------------------
    //6. generate info
    //info->Append(Form("Singal to Noise:   \t%d\n", id));
    info->Append(Form("Bary Center:       \t(%.2f, %.2f)\n", mBx, mBy));
    info->Append(Form("Principal Axis:  \tk=%.2f, b=%.2f\n", lPk, lPb));
    info->Append(Form("Conversion Point:  \t(%.2f, %.2f)\n", mCx, mCy));
    info->Append(Form("Ejection Axis:  \tk=%.2f, b=%.2f\n", lCk, lCb));
    info->Append(Form("Mom 3rd:        \t%.5f\n", mom3rd));
    info->Append(Form("Phi:            \t%.1f %.1f %.1f\n", k2, phimin,phimax));
}

//______________________________________________________________________________
void MyEventClass::Fill2DPlot(TH2F *h)
{
    for (int i = xmin; i <= xmax; i++)
        for (int j = ymin; j <= ymax; j++)
            h->SetBinContent(i + 1, j + 1, h->GetBinContent(i + 1, j + 1) + f2D->GetBinContent(i + 1, j + 1));
}

void MyEventClass::FillBaryCenter(TH2F *h)
{
    h->Fill(mBx, mBy);
}

void MyEventClass::FillIPpoint(TH2F *h)
{
    h->Fill(mCx, mCy);
}

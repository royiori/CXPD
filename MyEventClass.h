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
// fit function
//
double lineFunc(double *_x, double *par)
{
    double b = par[0];
    double k = par[1];
    double x = _x[0];
    double y = _x[1];
    double dist;

    dist = fabs(k * x + b - y) / sqrt(1 + k * k);
    return 1 / dist;
}

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
    double x, y;
    double dist;

    double k = p[1];
    double b = p[0];

    for (int i = 0; i < n; ++i)
    {
        if (values[i] == 0)
            continue;
        if (p[0] == 0 && p[1] == 0)
            continue;

        x = coords[i].first;
        y = coords[i].second;

        dist = fabs(k * x + b - y) / sqrt(1 + k * k);

        chi2 += values[i] * dist * dist;
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
        dataQFlag = 0;
        id = _id;
        xmin = _xmin;
        xmax = _xmax;
        ymin = _ymin;
        ymax = _ymax;
        nx = xmax - xmin + 1;
        ny = ymax - ymin + 1;
        f2D = NULL;
        mBcenter = NULL;
        lPrinAxis = NULL;
        mCovPoint = NULL;
        lCovAxis = NULL;
        info = NULL;

        data.resize((nx));
        for (int i = 0; i < nx; i++)
            data[i].resize(ny);
    }

    virtual ~MyEventClass();

    void Fill(int x, int y, double _b) { data[x - xmin][y - ymin] = _b; };
    void GenerateHist(TH2F *ped);
    void AnalysisHist(TH2F *ped);
    void Fill2DPlot(TH2F*);

    TH2F *Get2DPlot() { return f2D; };
    TH2F *Get2DRawPlot() { return f2D_raw; };

    TMarker *GetBaryCenterAsMarker() { return mBcenter; }
    TLine *GetPrincipalAxis() { return lPrinAxis; }
    TMarker *GetConvertionPoint() { return mCovPoint; }
    TLine *GetCovertionAxis() { return lCovAxis; }
    TString *GetInfo() { return info; };

    int GetDataQuality() {return dataQFlag;}
  private:
    int xmin;
    int xmax;
    int ymin;
    int ymax;
    int nx;
    int ny;
    int id;

    int clusterSize;
    int dataQFlag; //=0 means no data, =1 means has data
    double pulseHeight;

    double mBx, mBy;
    double lPk, lPb;
    double mCx, mCy;
    double lCk, lCb;

    TH2F *f2D;
    TH2F *f2D_raw;
    
    TMarker *mBcenter;
    TMarker *mCovPoint;
    TLine *lPrinAxis;
    TLine *lCovAxis;
    TString *info;

    vector<vector<double>> data;

    void EtchHistogram(TH2F *, TH2F *);
    void ExpandHistogram(TH2F *, TH2F *);
};

MyEventClass::~MyEventClass()
{
    // Destructor.
}

//______________________________________________________________________________
void MyEventClass::EtchHistogram(TH2F *f0, TH2F *f1)
{
    const int netch = 3;
    int ietch = 1;
    int etchData[netch][netch] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    //... etching
    for(int i=xmin; i<nx; i++)
        for(int j=ymin; j<ny; j++)
            f1->SetBinContent(i+1, j+1, 0);

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
    for(int i = xmin; i < nx; i++)
        for(int j = ymin; j < ny; j++)
            f1->SetBinContent(i+1, j+1, 0);

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
void MyEventClass::GenerateHist(TH2F *hPed)
{
    double xped = 0;
    double yped = 0;
    double siz = 1;

    if (f2D != NULL)
    {
        delete f2D;
        delete f2D_raw;
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
            f2D_raw->Fill((i - xped) * siz, (j - yped) * siz, (q < 0) ? 0 : q);

            //.. stored for fit..
            coords.push_back(std::make_pair((i - xped) * siz, (j - yped) * siz));
            values.push_back(q);
            errors.push_back(sqrt(q));
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
                f2D->SetBinContent(i+1, j+1, f2D_raw->GetBinContent(i+1, j+1));
            }

    delete ftmp0;
    delete ftmp1;

    //... generate info
    if (info != NULL)
        delete info;
    info = new TString();
    info->Append(Form("Event Number:      \t%d\n", id));
    info->Append(Form("Cluster Size:         \t%d\n", clusterSize));
    info->Append(Form("Pulse Height:      \t%.2f\n", pulseHeight));

    if (clusterSize < 20 || pulseHeight < 200)
        return;
    dataQFlag = 1;
}

void MyEventClass::AnalysisHist(TH2F *ped)
{
    //.. the barycenter point ..
    mBx = f2D->GetMean(1);
    mBy = f2D->GetMean(2);
    if (mBcenter != NULL)
        delete mBcenter;
    mBcenter = new TMarker(mBx, mBy, 30);
    mBcenter->SetMarkerColor(kRed);
    mBcenter->SetMarkerSize(2.6);

    //.. fit for principle line ..
    TF2 *func = new TF2("func", lineFunc, -1, 1, -1, 1, 2);

    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter *minuit = TVirtualFitter::Fitter(0, 10);
    minuit->SetParameter(0, func->GetParName(0), func->GetParameter(0), 0.01, 0, 0);
    minuit->SetParameter(1, func->GetParName(1), func->GetParameter(1), 0.01, 0, 0);
    minuit->SetFCN(myFcn);

    double arglist[100];
    arglist[0] = 0;
    minuit->ExecuteCommand("SET PRINT", arglist, 2);

    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.01; // tolerance
    minuit->ExecuteCommand("MIGRAD", arglist, 2);

    //get result
    double minParams[2];
    double parErrors[2];
    for (int i = 0; i < 2; ++i)
    {
        minParams[i] = minuit->GetParameter(i);
        parErrors[i] = minuit->GetParError(i);
    }

    double chi2, edm, errdef;
    int nvpar, nparx;
    minuit->GetStats(chi2, edm, errdef, nvpar, nparx);

    lPb = minParams[0];
    lPk = minParams[1];
    if (lPrinAxis != NULL)
        delete lPrinAxis;
    double xmin = f2D->GetXaxis()->GetXmin();
    double xmax = f2D->GetXaxis()->GetXmax();
    double ymin = f2D->GetYaxis()->GetXmin();
    double ymax = f2D->GetYaxis()->GetXmax();
    double _xmin = (ymin - lPb) / lPk;
    double _xmax = (ymax - lPb) / lPk;
    double _ymin = lPk * xmin + lPb;
    double _ymax = lPk * xmax + lPb;
    std::vector<std::pair<double, double>> _line;
    if (xmin <= _xmin && _xmin < xmax)
        _line.push_back(std::make_pair(_xmin, ymin));
    if (ymin <= _ymin && _ymin < ymax)
        _line.push_back(std::make_pair(xmin, _ymin));
    if (xmin <= _xmax && _xmax < xmax)
        _line.push_back(std::make_pair(_xmax, ymax));
    if (ymin <= _ymax && _ymax < ymax)
        _line.push_back(std::make_pair(xmax, _ymax));
    lPrinAxis = new TLine(_line[0].first, _line[0].second, _line[1].first, _line[1].second);
    lPrinAxis->SetLineColor(kRed);

    //.. fit for ejection line ..
    double max = 0;
    int imax = -1;

    double mom2nd = 0, mom3rd = 0;
    for (unsigned long i = 0; i < coords.size(); i++)
    {
        double x = coords[i].first;
        double y = coords[i].second;
        double q = values[i];
        if (q == 0)
            continue;

        double d2 = q * ((x - mBx) * (x - mBx) + (y - mBy) * (y - mBy));
        double d = (y - lPb - lPk * x) / sqrt(1 + lPk * lPk);
        double d3 = d * d * d;

        mom2nd += d2;
        mom3rd += d3;
        if (max < d2)
        {
            max = d2;
            imax = i;
        }
    }

    cout << "Max 2nd moment: (" << coords[imax].first << ", " << coords[imax].second << ")" << endl;

    if (mCovPoint != NULL)
        delete mCovPoint;
    mCx = coords[imax].first;
    mCy = coords[imax].second;
    mCovPoint = new TMarker(mCx, mCy, 30);
    mCovPoint->SetMarkerColor(kBlue);
    mCovPoint->SetMarkerSize(2);

    double sumx = 0;
    double sumy = 0;

    for (unsigned long i = 0; i < coords.size(); i++)
    {
        double xx = coords[i].first - mCx;
        double yy = coords[i].second - mCy;
        double dd = sqrt(xx * xx + yy * yy);
        double q = values[i];
        if (dd == 0 || q == 0)
            continue;
        //cout<<i<<" "<<xx<<" "<<yy<<" "<<dd<<"\n";
        //if (dd < 3 * siz || dd > 5 * siz)
        //    continue;

        double frac = q / dd; //weighted as Q/d
        sumx += frac * xx / dd;
        sumy += frac * yy / dd;
    }

    cout << "sum: " << sumy << " " << sumx << endl;

    lCk = sumy / sumx;
    lCb = coords[imax].second - lCk * coords[imax].first;

    cout << "y = " << lCk << "*x + " << lCb << endl;

    if (lCovAxis != NULL)
        delete lCovAxis;
    lCovAxis = new TLine((-1 - lCb) / lCk, -1, (1 - lCb) / lCk, 1);
    lCovAxis->SetLineColor(kBlue);

    //info->Append(Form("Singal to Noise:   \t%d\n", id));
    info->Append(Form("Bary Center:       \t(%.2f, %.2f)\n", mBx, mBy));
    info->Append(Form("Principal Axis:  \tk=%.2f, b=%.2f\n", lPk, lPb));
    info->Append(Form("Conversion Point:  \t(%.2f, %.2f)\n", mCx, mCy));
    info->Append(Form("Ejection Axis:  \tk=%.2f, b=%.2f\n", lCk, lCb));
}

//______________________________________________________________________________
void MyEventClass::Fill2DPlot(TH2F *h)
{
    for (int i = xmin; i <= xmax; i++)
        for (int j = ymin; j <= ymax; j++)
            h->SetBinContent(i+1, j+1, h->GetBinContent(i+1, j+1) + f2D->GetBinContent(i+1, j+1));
}

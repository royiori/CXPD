
#include "MyEventClass.h"

const double THRESHOLD = 0;

using namespace std;

//______________________________________________________________________________
// 用阶矩方法重建时所调用的拟合参数及函数
//----
// chisquare function
//
vector<pair<double, double>> coords;
vector<double> values; //Q
vector<double> errors; //dQ

int llm = 0;
double mBx = 0, mBy = 0, qtot = 0;    //bary
double mBx2 = 0, mBy2 = 0, qtot2 = 0; //impact
double mom3rd, qDistinguish;
double lPk, lPb;
double k2, b2;

void FcnToFitBaryLine(Int_t & /*nPar*/, Double_t * /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)
{
    int n = coords.size();
    double chi2 = 0;
    double x, y;
    double dist;

    double k = p[0];
    double b = mBy - k * mBx;
    //cout<<"hei "<<mBy<<" "<<mBx<<endl;
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

void FcnToFitEjectLine(Int_t & /*nPar*/, Double_t * /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)
{
    int n = coords.size();
    double chi2 = 0;
    double x, y;
    double dist1;
    double dist2;

    double k = p[0];
    double b = mBy2 - k * mBx2;

    for (int i = 0; i < n; ++i)
    {
        if (values[i] == 0)
            continue;

        x = coords[i].first;
        y = coords[i].second;
        double dNeed = (lPk == 0) ? (x - mBx) : (y - b2 - k2 * x) / sqrt(1 + k2 * k2); //去除一半的操作
        if (dNeed * mom3rd < 0)
            continue;

        dist1 = fabs(k * x + b - y) / sqrt(1 + k * k);
        dist2 = sqrt((x - mBx2) * (x - mBx2) + (y - mBy2) * (y - mBy2));
        if (dist2 != 0)
            chi2 += values[i] * dist1 * dist1;
    }
    fval = (chi2 == 0) ? 1E19 : chi2;
}

//______________________________________________________________________________
// 用Bessel曲线重建时所调用的函数， 三阶Bessel曲线
// 图像解释：https://www.cnblogs.com/jay-dong/archive/2012/09/26/2704188.html

double FcnOfBesselLine(Double_t t, vector<pair<double, double>> plist, int XY)
{
    if (t < 0 || t > 1 || (int)plist.size() != 4)
        return 0;

    double p0 = (XY == 0) ? plist[0].first : plist[0].second;
    double p1 = (XY == 0) ? plist[1].first : plist[1].second;
    double p2 = (XY == 0) ? plist[2].first : plist[2].second;
    double p3 = (XY == 0) ? plist[3].first : plist[3].second;

    return p0 * pow(1 - t, 3) + 3 * p1 * t * pow(1 - t, 2) + 3 * p2 * t * t * (1 - t) + p3 * pow(t, 3);
}

double BX(Double_t t, vector<pair<double, double>> plist) { return FcnOfBesselLine(t, plist, 0); }
double BY(Double_t t, vector<pair<double, double>> plist) { return FcnOfBesselLine(t, plist, 1); }

//画bessel曲线
void DrawBesselLine(vector<pair<double, double>> plist, double ts = 0.0001, const char *opt = "")
{
    for (double t = 0; t < 1 - ts; t += ts)
    {
        TLine *l = new TLine(BX(t, plist), BY(t, plist), BX(t + ts, plist), BY(t + ts, plist));
        l->SetLineColor(kRed);
        l->Draw(opt);
    }
}

void FcnToFitBaryLineByBesselLine(Int_t & /*nPar*/, Double_t * /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)
{
    int n = coords.size();
    double chi2 = 0;
    double x, y;

    //(p[0], p[1]) & (p[6], p[7])必须在cluster内部

    //
    vector<pair<double, double>> plist;
    plist.push_back(make_pair(p[0], p[1]));
    plist.push_back(make_pair(p[2], p[3]));
    plist.push_back(make_pair(p[4], p[5]));
    plist.push_back(make_pair(p[6], p[7]));
    double ts = p[8]; //bessel的扫描精度，请固定住不要变化

    //cout<<"hei "<<mBy<<" "<<mBx<<endl;
    vector<double> bxlist;
    vector<double> bylist;
    for (double t = 0; t < 1 - ts; t += ts)
    {
        double bx = BX(t, plist);
        double by = BY(t, plist);
        bxlist.push_back(bx);
        bylist.push_back(by);
    }

    for (int i = 0; i < n; ++i)
    {
        if (values[i] == 0)
            continue;

        x = coords[i].first;
        y = coords[i].second;
        double distMin = 9999;

        for (int j = 0; j < (int)bxlist.size(); j++)
        {
            double bx = bxlist[j];
            double by = bylist[j];
            double dtmp = sqrt((x - bx) * (x - bx) + (y - by) * (y - by));
            distMin = (dtmp < distMin) ? dtmp : distMin;
        }

        chi2 += values[i] * distMin * distMin;
    }

    fval = (chi2 == 0) ? 1E19 : chi2;
}

//______________________________________________________________________________
// 构造/析构函数
MyEventClass::MyEventClass(int _id, int _xmin, int _xmax, int _ymin, int _ymax)
{
    dataQFlag = QF_PED;
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

    defVal = -1000;
    mBx = mBy = lPk = lPb = defVal;
    mCx = mCy = lCk = lCb = defVal;
    aTheta1 = aTheta2 = defVal;

    nEtchingMatrix = 3;
    nExpandMatrix = 3;
    nEtchExpand = 2;
    nFattyX = 6;
    nFattyY = 6;
    nEllipticity = 1.5;
    rMinScale = 1.;
    rMaxScale = 1.;

    data.resize((nx));
    for (int i = 0; i < nx; i++)
        data[i].resize(ny);
}

MyEventClass::~MyEventClass()
{
    // Destructor.
    if (f2D != NULL)
        delete f2D;
    if (f2D_raw != NULL)
        delete f2D_raw;
    if (mBcenter != NULL)
        delete mBcenter;
    if (mCovPoint != NULL)
        delete mCovPoint;
    if (lPrinAxis1 != NULL)
        delete lPrinAxis1;
    if (lPrinAxis2 != NULL)
        delete lPrinAxis2;
    if (e1 != NULL)
        delete e1;
    if (e2 != NULL)
        delete e2;
    if (lCovAxis != NULL)
        delete lCovAxis;
    if (info != NULL)
        delete info;

    data.clear();
}

//______________________________________________________________________________
// 腐蚀/膨胀处理
void MyEventClass::EtchHistogram(TH2F *f0, TH2F *f1)
{
    int netch = nEtchingMatrix;
    int ietch = netch / 2.;
    netch = (netch < 2) ? 2 : netch;
    netch = (netch > 5) ? 5 : netch;

    vector<vector<int>> etchData;
    etchData.resize(netch);
    for (int i = 0; i < netch; i++)
    {
        etchData[i].resize(netch);
        for (int j = 0; j < netch; j++)
            etchData[i][j] = 1;
    }

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
                    if (m >= netch || n >= netch || etchData[m][n] == 0)
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
    int netch = nExpandMatrix;
    int ietch = netch / 2.;
    netch = (netch < 2) ? 2 : netch;
    netch = (netch > 5) ? 5 : netch;

    vector<vector<int>> etchData;
    etchData.resize(netch);
    for (int i = 0; i < netch; i++)
    {
        etchData[i].resize(netch);
        for (int j = 0; j < netch; j++)
            etchData[i][j] = 1;
    }

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
                    if (m >= netch || n >= netch || etchData[m][n] == 0)
                        continue;
                    f1->SetBinContent(ii + 1, jj + 1, 1);
                }
        }
    return;
}

//______________________________________________________________________________
// 0. 通用的处理函数
void MyEventClass::GenerateHist(TH2F *hPed, bool anaflag)
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
    f2D->GetXaxis()->SetTitle("X Pixel(every pixel=0.08mm)");
    f2D->GetYaxis()->SetTitle("Y pixel(every pixel=0.08mm)");

    //... store to data
    coords.clear();
    values.clear();
    errors.clear();
    for (int i = xmin; i < nx; i++)
        for (int j = ymin; j < ny; j++)
        {
            double ped = THRESHOLD;
            if (hPed != NULL)
            {
                ped = hPed->GetBinContent(i + 1, j + 1);
            }

            double q = data[i - xmin][j - ymin] - ped;
            q = (q < 0) ? 0 : q;
            f2D_raw->Fill(i + 1, j + 1, q);
        }

    //... deal histogram
    TH2F *ftmp0 = new TH2F("ftmp0", "ftmp", nx, 1, xmax - xmin + 1, ny, 1, ymax - ymin + 1);
    TH2F *ftmp1 = new TH2F("ftmp1", "ftmp", nx, 1, xmax - xmin + 1, ny, 1, ymax - ymin + 1);

    if (nEtchExpand > 0)
    {
        EtchHistogram(f2D_raw, ftmp0);
        ExpandHistogram(ftmp0, ftmp1);

        for (int i = 1; i < nEtchExpand; i++)
        {
            ExpandHistogram(ftmp1, ftmp0);
            EtchHistogram(ftmp0, ftmp1);
        }
    }

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
                coords.push_back(make_pair(i + 1, j + 1));
                values.push_back(f2D_raw->GetBinContent(i + 1, j + 1));
                errors.push_back(sqrt(f2D_raw->GetBinContent(i + 1, j + 1)));
            }

    delete ftmp0;
    delete ftmp1;

    //... generate info
    info = new TString();
    info->Append(Form("Event Number:      \t%d\n", id));
    info->Append(Form("Cluster Size:      \t%d\n", clusterSize));
    info->Append(Form("Pulse Height:      \t%.2f\n", pulseHeight));

    // No signal contained, treated as a pedstral
    if (clusterSize < 20 || pulseHeight < 200)
    {
        dataQFlag = QF_PED;
        f2D->SetTitle(Form("Event %d (identified as Ped)", id));
        return;
    }

    // Signal found, need to do a further analysis
    dataQFlag = QF_HIT;

    if (anaflag)
        AnalysisHist2();

    return;
}

void MyEventClass::Fill2DPlot(TH2F *h)
{
    for (int i = xmin; i <= xmax; i++)
        for (int j = ymin; j <= ymax; j++)
            h->SetBinContent(i + 1, j + 1, h->GetBinContent(i + 1, j + 1) + f2D->GetBinContent(i + 1, j + 1));
}

//______________________________________________________________________________
// 1. 算法1
//    采用阶矩来进行重建的算法
void MyEventClass::AnalysisHist()
{
    //----------------------
    //0. init

    int N = coords.size();

    double arglist[100];
    double chi2, edm, errdef;
    int nvpar, nparx;

    //----------------------
    //1. find barycenter
    {
        mBx = 0;
        mBy = 0;
        qtot = 0;

        for (int i = 0; i < N; i++)
        {
            double x = coords[i].first;
            double y = coords[i].second;
            double q = values[i];

            //cout<<i<<": "<<x<<", "<<y<<"  - "<<q<<endl;

            qtot += q;
            mBx += q * x;
            mBy += q * y;
        }

        mBx /= qtot;
        mBy /= qtot;

        //1.1 draw result
        mBcenter = new TMarker(mBx, mBy, 30);
        mBcenter->SetMarkerColor(kRed);
        mBcenter->SetMarkerSize(2.6);
        mBcenter->Draw();
    }

    //----------------------
    //2. fit barycenter-line
    {
        TVirtualFitter::SetDefaultFitter("Minuit");
        TVirtualFitter *minuit = TVirtualFitter::Fitter(0, 10);
        minuit->SetParameter(0, "k", 0, 0.01, 0, 0); //0.01
        minuit->SetFCN(FcnToFitBaryLine);

        arglist[0] = 0;
        minuit->ExecuteCommand("SET NOWarnings", arglist, 1);

        arglist[0] = -1;
        minuit->ExecuteCommand("SET PRINT", arglist, 1);

        arglist[0] = 5000;  // number of function calls
        arglist[1] = 0.001; // tolerance
        minuit->ExecuteCommand("MIGRAD", arglist, 2);

        //2.1 get result
        minuit->GetStats(chi2, edm, errdef, nvpar, nparx);

        lPk = minuit->GetParameter(0);
        lPb = mBy - lPk * mBx;
        if (abs(lPk) < 0.01)
            lPk = 0;

        //2.2 draw result
        double lx1 = (lPk == 0) ? xmin : (ymin - lPb) / lPk;
        double ly1 = (lPk == 0) ? mBy : ymin;
        double lx2 = (lPk == 0) ? xmax : (ymax - lPb) / lPk;
        double ly2 = (lPk == 0) ? mBy : ymax;
        aTheta1 = atan(lPk);

        lPrinAxis1 = new TLine(lx1, ly1, lx2, ly2);
        lPrinAxis1->SetLineColor(kBlack);
        lPrinAxis1->Draw();
    }

    //----------------------
    //2.5 cut the fatty events
    {
        double minx = 99, miny = 99;
        double maxx = -1, maxy = -1;

        for (int i = 0; i < N; i++)
        {
            double x = coords[i].first;
            double y = coords[i].second;

            minx = (minx > x) ? x : minx;
            miny = (miny > y) ? y : miny;
            maxx = (maxx < x) ? x : maxx;
            maxy = (maxy < y) ? y : maxy;
        }

        if (maxx - minx < nFattyX && maxy - miny < nFattyY)
        {
            if ((maxx - minx) / (maxy - miny) < nEllipticity && (maxx - minx) / (maxy - miny) > 1 / nEllipticity)
            { //判断fat事例的圆形程度
                dataQFlag = QF_FAT;
                f2D->SetTitle(Form("Event %d (identified as Fatty events)", id));

                mCovPoint = new TMarker(mBx, mBy, 31);
                mCovPoint->SetMarkerColor(kRed);
                mCovPoint->SetMarkerSize(2.6);
                mCovPoint->Draw();

                lCovAxis = new TLine(lPrinAxis1->GetX1(), lPrinAxis1->GetY1(), lPrinAxis1->GetX2(), lPrinAxis1->GetY2());
                lCovAxis->SetLineColor(kRed);
                lCovAxis->Draw();
                return;
            }
        }
    }

    //----------------------
    //3. cal mom3rd
    {
        k2 = (lPk == 0) ? 1 : -1 / lPk;
        b2 = (lPk == 0) ? 1 : mBy - k2 * mBx;

        //3.1 find suitable point for mom3rd
        mom3rd = 0;
        for (int i = 0; i < N; i++)
        {
            double x = coords[i].first;
            double y = coords[i].second;
            double q = values[i];
            double d = (lPk == 0) ? (x - mBx) : (y - k2 * x - b2) / sqrt(1 + k2 * k2);
            double d3 = d * d * d * q;
            mom3rd += d3;
        }
    }

    //----------------------
    //4. calculate 2nd mom max and define rmin
    double theta2 = (lPk == 0) ? TMath::Pi() / 2 : atan(-1 / lPk);

    double mom2ndMid0 = 0;
    double mom2ndMid1 = 0;
    double mom2ndMid2 = 0;
    double shapeRatio = 0;
    double rmin;
    double rmax = 0;

    int iterations = 0;
    double chooseIteration[4];
    double intermediateResult = 0;
    double caluDecent = 0;
    double pureMom2nd;
    double calumom2nd;
    double stepB2 = (mom3rd < 0) ? -0.08 / cos(theta2) : 0.08 / cos(theta2);
    double clusterlength = 0;

    const int iteranum = 10;

    do
    {
        mom2nd = 0;
        mom2ndMax = 0;
        mom2ndMin = 0;
        rmin = 0;
        pureMom2nd = 0;
        for (int i = 0; i < N; i++)
        {
            double x = coords[i].first;
            double y = coords[i].second;
            double q = values[i];

            double x0 = (lPk * lPk * mBx + lPk * (y - mBy) + x) / (lPk * lPk + 1); //到主轴的mom2nd
            double y0 = lPk * (x0 - mBx) + mBy;

            double d = (lPk == 0) ? (x - mBx) : (y - b2 - k2 * x) / sqrt(1 + k2 * k2);
            if (d * mom3rd < 0)
                continue;

            mom2ndMid0 = (pow(x - x0, 2) + pow(y - y0, 2));
            mom2ndMid1 = q * (pow(mBx - x0, 2) + pow(mBy - y0, 2));
            mom2ndMid2 = q * mom2ndMid0;

            mom2nd += mom2ndMid1 + mom2ndMid2;
            mom2ndMax += mom2ndMid1;
            mom2ndMin += mom2ndMid2;
            pureMom2nd += mom2ndMid2;

            rmin += (mom2ndMid2);
        }
        if (iterations == 0)
        {
            calumom2nd = pureMom2nd;
            clusterlength = mom2ndMax / mom2ndMin;
        }
        b2 = (clusterlength > 5) ? b2 + stepB2 : b2 - stepB2;
        caluDecent = abs(pureMom2nd - intermediateResult) / calumom2nd;
        //cout<<"the "<<iterations<<" times "<<pureMom2nd<<" value "<<caluDecent<<endl;
        if (iterations == 0 || (caluDecent < 0.05 && chooseIteration[0] == 0 && pureMom2nd > 20))
        { //限制条件找合适的截取矩
            chooseIteration[0] = iterations;
            chooseIteration[1] = mom2ndMax;
            chooseIteration[2] = mom2ndMin;
            chooseIteration[3] = rmin;
        }
        intermediateResult = pureMom2nd;
        iterations++;
    } while (iterations < iteranum);

    b2 = (clusterlength > 5) ? b2 - (iteranum - chooseIteration[0]) * stepB2 : b2 + (iteranum - chooseIteration[0]) * stepB2;
    mom2ndMax = chooseIteration[1];
    mom2ndMin = chooseIteration[2];
    rmin = chooseIteration[3];

    double xintersecton = (b2 - lPb) / (lPk - k2);
    double yintersecton = k2 * xintersecton + b2;

    //4.2 draw result
    double lx1 = (lPk == 0) ? mBx : (ymin - b2) / k2;
    double ly1 = ymin;
    double lx2 = (lPk == 0) ? mBx : (ymax - b2) / k2;
    double ly2 = ymax;
    lPrinAxis2 = new TLine(lx1, ly1, lx2, ly2);
    lPrinAxis2->SetLineColor(kBlue);
    lPrinAxis2->Draw();

    //mom2nd /=  qtot;
    mom2ndMax /= qtot;
    mom2ndMin /= qtot;
    shapeRatio = (mom2ndMin == 0) ? 0 : mom2ndMax / mom2ndMin;

    rmin /= 2 * qtot;
    rmin = sqrt(rmin);
    rmin *= (rMinScale > 0) ? rMinScale : 1.;
    rmax = rmin * 9;
    rmax *= (rMaxScale > 0) ? rMaxScale : 1.;

    //4.1 draw TEllipse
    double phimin = theta2 / TMath::Pi() * 180;
    double phimax = theta2 / TMath::Pi() * 180;
    if (mom3rd < 0)
        phimin -= 180;
    if (mom3rd > 0)
        phimax += 180;
    //    if (k2 < 0 && mom3rd > 0)
    //	    phimin -= 180;
    //    if (k2 < 0 && mom3rd < 0)
    //        phimax += 180;

    e1 = new TEllipse(xintersecton, yintersecton, rmin, rmin, phimin, phimax);
    e2 = new TEllipse(xintersecton, yintersecton, rmax, rmax, phimin, phimax);
    e1->SetFillColorAlpha(kWhite, 0);
    e1->SetFillStyle(0);
    e2->SetFillColorAlpha(kWhite, 0);
    e2->SetFillStyle(0);
    e1->Draw();
    e2->Draw();

    //----------------------
    //5. find IP
    mBx2 = 0;
    mBy2 = 0;
    qtot2 = 0;

    for (int i = 0; i < N; i++)
    {
        double x = coords[i].first;
        double y = coords[i].second;
        double q = values[i];

        double d = (lPk == 0) ? (x - mBx) : (y - b2 - k2 * x) / sqrt(1 + k2 * k2);

        if (d * mom3rd < 0)
            continue;
        if (sqrt(pow(xintersecton - x, 2) + pow(yintersecton - y, 2)) < rmin)
            continue;
        if (sqrt(pow(xintersecton - x, 2) + pow(yintersecton - y, 2)) > rmax)
            continue;

        mBx2 += q * x;
        mBy2 += q * y;
        qtot2 += q;
    }

    mBx2 /= qtot2;
    mBy2 /= qtot2;
    mCx = mBx2;
    mCy = mBy2;

    //5.1 draw result
    mCovPoint = new TMarker(mCx, mCy, 31);
    mCovPoint->SetMarkerColor(kRed);
    mCovPoint->SetMarkerSize(2.6);
    mCovPoint->Draw();

    //----------------------
    //6. Fit eject line
    TVirtualFitter *minuit2 = TVirtualFitter::Fitter(0, 10);
    minuit2->SetParameter(0, "k", 0, 0.01, 0, 0);
    minuit2->SetFCN(FcnToFitEjectLine);

    arglist[0] = 0;
    minuit2->ExecuteCommand("SET NOWarnings", arglist, 1);

    arglist[0] = -1;
    minuit2->ExecuteCommand("SET PRINT", arglist, 1);

    arglist[0] = 5000;  // number of function calls
    arglist[1] = 0.001; // tolerance
    minuit2->ExecuteCommand("MIGRAD", arglist, 2);

    //5.1 get result
    minuit2->GetStats(chi2, edm, errdef, nvpar, nparx);

    double k3 = minuit2->GetParameter(0);
    double b3 = mBy2 - k3 * mBx2;
    lCk = k3;
    lCb = b3;
    aTheta2 = atan(lCk);

    //5.2 draw result
    lx1 = (k3 == 0) ? mBx2 : (ymin - b3) / k3;
    ly1 = ymin;
    lx2 = (k3 == 0) ? mBx2 : (ymax - b3) / k3;
    ly2 = ymax;

    lCovAxis = new TLine(lx1, ly1, lx2, ly2);
    lCovAxis->SetLineColor(kRed);
    lCovAxis->Draw();

    if (aTheta1 > 0)
        aTheta1 = (mom3rd > 0) ? aTheta1 - TMath::Pi() : aTheta1;
    if (aTheta1 < 0)
        aTheta1 = (mom3rd < 0) ? aTheta1 + TMath::Pi() : aTheta1;

    if (aTheta2 > 0)
        aTheta2 = (mom3rd > 0) ? aTheta2 - TMath::Pi() : aTheta2;
    if (aTheta2 < 0)
        aTheta2 = (mom3rd < 0) ? aTheta2 + TMath::Pi() : aTheta2;
    //k3 +=0;
    if (mom3rd == 0 || dataQFlag == 2)
    {
        aTheta1 = 2 * TMath::Pi() - 0.05;
        aTheta2 = 2 * TMath::Pi() - 0.05;
        llm++;
    }

    //-------------------
    //6. generate info
    //info->Append(Form("Singal to Noise:   \t%d\n", id));
    info->Append(Form("Bary Center:       \t(%.2f, %.2f)\n", mBx, mBy));
    info->Append(Form("Principal Axis:  \tk=%.2f, b=%.2f\n", lPk, lPb));
    info->Append(Form("Conversion Point:  \t(%.2f, %.2f)\n", mCx, mCy));
    info->Append(Form("Ejection Axis:  \tk=%.2f, b=%.2f\n", lCk, lCb));
    info->Append(Form("Mom 2rd:        \t%.5f\n", mom2ndMax / mom2ndMin));
    info->Append(Form("Mom 3rd:        \t%.10f\n", mom3rd));
    info->Append(Form("r:        \t%.2f, \t%.2f\n", rmin, rmax));
    //info->Append(Form("aTheta1:   \t%.2f\n", aTheta1));
}

//_____________________________
// 采用阶矩来进行重建的中间变量填图
void MyEventClass::FillBaryCenter(TH2F *h)
{
    if (mBx != defVal && mBy != defVal)
        h->Fill(mBx, mBy);
}

void MyEventClass::FillIPpoint(TH2F *h)
{
    if (mCx != defVal && mCy != defVal)
        h->Fill(mCx, mCy);
}

void MyEventClass::FillPolarization1(TH1F *h)
{
    if (aTheta1 != defVal)
        h->Fill(aTheta1);
}

void MyEventClass::FillPolarization2(TH1F *h)
{
    if (aTheta2 != defVal)
        h->Fill(aTheta2);
}

void MyEventClass::Filllengththelongest(TH1F *h)
{
    h->Fill(sqrt(mom2nd));
}

//_____________________________
void MyEventClass::Draw2DResultMethod1(const char *opt)
{
    if (f2D != NULL)
        f2D->Draw(opt);
    if (mBcenter != NULL)
        mBcenter->Draw();
    if (mCovPoint != NULL)
        mCovPoint->Draw();
    if (lPrinAxis1 != NULL)
        lPrinAxis1->Draw();
    if (lCovAxis != NULL)
        lCovAxis->Draw();
}

void MyEventClass::DrawSearchResultMethod1()
{
    if (lPrinAxis2 != NULL)
        lPrinAxis2->Draw();
    if (e1 != NULL)
        e1->Draw();
    if (e2 != NULL)
        e2->Draw();
}

//______________________________________________________________________________
// 2. 算法2
//    采用bessel曲线描述径迹的算法
void MyEventClass::AnalysisHist2()
{
    int N = coords.size();

    double arglist[100];
    double chi2, edm, errdef;
    int nvpar, nparx;
    int imax = -1;
    int imin = -1;
    double xmin, xmax, ymin, ymax;

    plist.clear();

    //1. cluster分组寻找
    {
        double max = -1;
        for (int i = 0; i < N; i++)
            if (values[i] > max)
            {
                imax = i;
                max = values[i];
            }

        if (imax == -1)
            return;

        double dist = -1;
        for (int i = 0; i < N; i++)
            if (dist < sqrt(pow(coords[imax].first - coords[i].first, 2) + pow(coords[imax].second - coords[i].second, 2)))
            {
                dist = sqrt(pow(coords[imax].first - coords[i].first, 2) + pow(coords[imax].second - coords[i].second, 2));
                imin = i;
            }

        if (imin == -1)
            return;
    }

    //2. 拟合bessel曲线的形式
    {
        TVirtualFitter::SetDefaultFitter("Minuit");
        TVirtualFitter *minuit = TVirtualFitter::Fitter(0);
        minuit->SetParameter(0, "p0x", coords[imax].first, 0.01, 0, 0);
        minuit->SetParameter(1, "p0y", coords[imax].second, 0.01, 0, 0);
        minuit->SetParameter(2, "p1x", coords[imax].first, 0.01, 0, 0);
        minuit->SetParameter(3, "p1y", coords[imax].second, 0.01, 0, 0);
        minuit->SetParameter(4, "p2x", coords[imin].first, 0.01, 0, 0);
        minuit->SetParameter(5, "p2y", coords[imin].second, 0.01, 0, 0);
        minuit->SetParameter(6, "p3x", coords[imin].first, 0.01, 0, 0);
        minuit->SetParameter(7, "p3y", coords[imin].second, 0.01, 0, 0);
        minuit->SetParameter(8, "precision", 0.01, 0.01, 0, 0);
        minuit->FixParameter(8);
        minuit->SetFCN(FcnToFitBaryLineByBesselLine);

        arglist[0] = 0;
        minuit->ExecuteCommand("SET NOWarnings", arglist, 1);

        arglist[0] = -1;
        minuit->ExecuteCommand("SET PRINT", arglist, 1);

        arglist[0] = 5000;  // number of function calls
        arglist[1] = 0.001; // tolerance
        minuit->ExecuteCommand("MIGRAD", arglist, 2);
        //2.1 get result
        minuit->GetStats(chi2, edm, errdef, nvpar, nparx);

        plist.push_back(make_pair(minuit->GetParameter(0), minuit->GetParameter(1)));
        plist.push_back(make_pair(minuit->GetParameter(2), minuit->GetParameter(3)));
        plist.push_back(make_pair(minuit->GetParameter(4), minuit->GetParameter(5)));
        plist.push_back(make_pair(minuit->GetParameter(6), minuit->GetParameter(7)));

        //2.2 draw result
    }
}

void MyEventClass::Draw2DResultMethod2(const char *opt)
{
    if (plist.size() == 0)
        return;

    if (f2D != NULL)
        f2D->Draw(opt);
    DrawBesselLine(plist, 0.0001, opt);
}

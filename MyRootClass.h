#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"

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
   
   dist = fabs(k*x + b - y)/sqrt(1+k*k);
   return 1/dist;
}

//----
// chisquare function
//
std::vector<std::pair<double, double> > coords;
std::vector<double > values;
std::vector<double > errors;

void myFcn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
  int n = coords.size();
  double chi2 = 0;
  double x,y;
  double dist;

  double k = p[1];
  double b = p[0];

  for (int i = 0; i<n; ++i ) {
    if(values[i]==0) continue;
    if(p[0]==0 && p[1]==0) continue;

    x = coords[i].first;
    y = coords[i].second;

    dist = fabs(k*x + b - y)/sqrt(1+k*k);

    chi2 += values[i]*dist*dist;
  }
  fval = (chi2==0)?1E19:chi2;
}



//-----------------------------------------------------------------------------------------------
//
class MyEventClass
{
  public:
    MyEventClass(int _id, int _xmin, int _xmax, int _ymin, int _ymax)
    {
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

    virtual ~MyEventClass(){};

    void Fill(int x, int y, double _b) { data[x - xmin][y - ymin] = _b; };
    void GenerateHist(double ped);
    void FillPedestal(TH1F *h);

    TH2F *Get2DPlot() { return f2D; };
    TMarker *GetBaryCenterAsMarker() { return mBcenter; }
    TLine *GetPrincipalAxis() { return lPrinAxis; }
    TMarker *GetConvertionPoint() { return mCovPoint; }
    TLine *GetCovertionAxis() { return lCovAxis; }
    TString *GetInfo() {return info;};

    //private:
    int xmin;
    int xmax;
    int ymin;
    int ymax;
    int nx;
    int ny;
    int id;

    double mBx, mBy;
    double lPk, lPb;
    double mCx, mCy;
    double lCk, lCb;

    TH2F *f2D;
    TMarker *mBcenter;
    TMarker *mCovPoint;
    TLine   *lPrinAxis;
    TLine   *lCovAxis;
    TString *info;

    vector<vector<double>> data;
};

void MyEventClass::GenerateHist(double ped)
{
    int xped = 176;
    int yped = 150;
    double siz = 0.043;

    if(f2D!=NULL) delete f2D;
    f2D = new TH2F(Form("f2D_%d", id), "2D Plot", nx, (xmin-xped)*siz, (xmax-xped)*siz, ny, (ymin-yped)*siz, (ymax-yped)*siz);

    int cluster_size = 0;
    double pulse_height = 0;

    coords.clear();
    values.clear();
    errors.clear();
    for (int i = xmin; i <= xmax; i++)
    {
        for (int j = ymin; j <= ymax; j++)
        {
            double q = (data[i - xmin][j - ymin] - ped) < 0 ? 0 : data[i - xmin][j - ymin];
            f2D->Fill((i-xped)*siz, (j-yped)*siz, q);
            if(q>0) cluster_size++;
            pulse_height+=q;

            //.. stored for fit..
            coords.push_back( std::make_pair((i-xped)*siz, (j-yped)*siz));
            values.push_back( q );
            errors.push_back( sqrt(q));
        }
    }
    //.. the barycenter point ..
    mBx = f2D->GetMean(1);
    mBy = f2D->GetMean(2);
    if(mBcenter!=NULL) delete mBcenter;
    mBcenter = new TMarker(mBx, mBy, 30);
    mBcenter->SetMarkerColor(kRed);
    mBcenter->SetMarkerSize(2.6);

    //.. fit for principle line ..
    TF2 *func = new TF2("func",lineFunc,-1, 1,-1, 1, 2);

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
    if(lPrinAxis!=NULL) delete lPrinAxis;
    double xmin = f2D->GetXaxis()->GetXmin();
    double xmax = f2D->GetXaxis()->GetXmax();
    double ymin = f2D->GetYaxis()->GetXmin();
    double ymax = f2D->GetYaxis()->GetXmax();
    double _xmin = (ymin-lPb)/lPk;
    double _xmax = (ymax-lPb)/lPk;
    double _ymin = lPk*xmin+lPb;
    double _ymax = lPk*xmax+lPb;
    std::vector<std::pair<double, double> > _line;
    if(xmin <= _xmin && _xmin < xmax) _line.push_back( std::make_pair(_xmin, ymin)); 
    if(ymin <= _ymin && _ymin < ymax) _line.push_back( std::make_pair(xmin, _ymin));
    if(xmin <= _xmax && _xmax < xmax) _line.push_back( std::make_pair(_xmax, ymax)); 
    if(ymin <= _ymax && _ymax < ymax) _line.push_back( std::make_pair(xmax, _ymax));
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
        if(q==0) continue;

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

    if(mCovPoint!=NULL) delete mCovPoint;
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
        double  q = values[i];
        if(dd==0 || q==0) continue;
        //cout<<i<<" "<<xx<<" "<<yy<<" "<<dd<<"\n";
        //if (dd < 3 * siz || dd > 5 * siz)
        //    continue;

        double frac = q / dd; //weighted as Q/d
        sumx += frac * xx / dd;
        sumy += frac * yy / dd;
    }

    cout<<"sum: "<<sumy<<" "<<sumx<<endl;

    lCk = sumy / sumx;
    lCb = coords[imax].second - lCk * coords[imax].first;

    cout << "y = " << lCk << "*x + " << lCb << endl;

    if(lCovAxis!=NULL) delete lCovAxis;
    lCovAxis = new TLine((-1 - lCb) / lCk, -1, (1 - lCb) / lCk, 1);
    lCovAxis->SetLineColor(kBlue);

    if(info!=NULL) delete info;
    info = new TString();
    info->Append(Form("Event Number:      \t%d\n", id));
    info->Append(Form("Cluster Size:         \t%d\n", cluster_size));
    info->Append(Form("Pulse Height:      \t%.2f\n", pulse_height));
    //info->Append(Form("Singal to Noise:   \t%d\n", id));
    info->Append(Form("Bary Center:       \t(%.2f, %.2f)\n", mBx, mBy));
    info->Append(Form("Principal Axis:  \tk=%.2f, b=%.2f\n", lPk, lPb));
    info->Append(Form("Conversion Point:  \t(%.2f, %.2f)\n", mCx, mCy));
    info->Append(Form("Ejection Axis:  \tk=%.2f, b=%.2f\n", lCk, lCb));
}

void MyEventClass::FillPedestal(TH1F *h)
{
    for (int i = xmin; i <= xmax; i++)
    {
        for (int j = ymin; j <= ymax; j++)
        {
            h->Fill(data[i - xmin][j - ymin]);
        }
    }
}


//-----------------------------------------------------------------------------------------------
//
class MyRootClass
{
  public:
    MyRootClass()
    {
        ip = 0;
        ped = 0;

    };
    ~MyRootClass(){};

    void Analysis(TString);

    TString* GetInfo();
    TH2F *DrawPre();
    TH2F *DrawNext();
    TH1F *ButtonFunc1();
    void ButtonFunc2(double _ped) { ped = _ped; };
    TMarker *ButtonFunc3();
    TLine *ButtonFunc4();
    TMarker *ButtonFunc5();
    TLine *ButtonFunc6();

  private:
    int ip;
    double ped;
    TH1F *hPed;
    vector<MyEventClass *> fEventList;
};

void MyRootClass::Analysis(TString filePath)
{
    ip = 0;
    ped = 0;
    fEventList.clear();

    cout << "--> Opening: " << filePath << endl;
    ifstream ifSignal(filePath, ios::binary);

    int nEvent = 0;
    short xmin, xmax, ymin, ymax;
    unsigned short hdata[10];
    int nByte = 0;

    while (ifSignal.good())
    {
        unsigned short _data;
        ifSignal.read((char *)(&_data), sizeof(_data));
        hdata[0] = _data;
        nByte++;

        if (_data != 65535)
        {
            cout << "bad: " << _data << " at " << nByte << endl;
            continue;
        }
        else
        {
            for (int j = 1; j <= 9; j++)
            {
                ifSignal.read((char *)(&_data), sizeof(_data));
                hdata[j] = _data;
                nByte++;
            }
            ymin = hdata[1];
            ymax = hdata[2];
            xmin = hdata[3];
            xmax = hdata[4];
            //if (xmin>300||xmax>300) continue;
            //if (ymin>352||ymax>352) continue;
            if (nEvent > 12590)
                cout << "----> x: " << xmin << " " << xmax << " y: " << ymin << " " << ymax << endl;
            if (nEvent > 12590)
                cout << "----> evtid :" << hdata[5] << endl;
            if (nEvent > 12590)
                cout << "----> timestamp: " << hdata[6] << " " << hdata[7] << " " << hdata[8] << " " << hdata[9] << endl;

            MyEventClass *fEvent = new MyEventClass(nEvent, xmin, xmax, ymin, ymax);

            for (int i = xmin; i <= xmax; i++)
            {
                for (int j = ymin; j <= ymax; j++)
                {
                    ifSignal.read((char *)(&_data), sizeof(_data));
                    fEvent->Fill(i, j, _data);
                    nByte++;
                }
            }
            fEvent->GenerateHist(ped);
            fEventList.push_back(fEvent);
        }
        nEvent++;
        if (nEvent % 100 == 0)
            cout << "Dealed: " << nEvent << endl;
        if (nEvent > 100)
            break;
    }

    cout << "----> " << (unsigned long)ifSignal.tellg() << " is readed." << endl;
    cout << "----> " << nEvent << " events recorded." << endl;
    ifSignal.close();
}

TH2F *MyRootClass::DrawPre()
{
    if (ip <= 0)
        return NULL;
    ip--;
    fEventList.at(ip)->GenerateHist(ped);
    return fEventList.at(ip)->Get2DPlot();
}

TH2F *MyRootClass::DrawNext()
{
    if ((unsigned long)ip == fEventList.size() - 1)
        return NULL;
    ip++;
    fEventList.at(ip)->GenerateHist(ped);
    return fEventList.at(ip)->Get2DPlot();
}

TString *MyRootClass::GetInfo()
{
    return fEventList.at(ip)->GetInfo();
}

TH1F *MyRootClass::ButtonFunc1()
{
    TString hname("hPed");
    if (gDirectory->Get(hname) != NULL)
    {
        hPed = (TH1F *)gDirectory->Get(hname);
        delete hPed;
    }
    hPed = new TH1F(hname, "Pedestal", 100, 0, 100);
    for (unsigned long ip = 0; ip < fEventList.size(); ip++)
    {
        fEventList.at(ip)->FillPedestal(hPed);
    }
    TF1 *fS = new TF1("fS", "gaus", 0, 180);
    fS->SetParameter(1, 0);
    fS->SetParameter(2, hPed->GetRMS());
    hPed->Fit(fS, "RB");
    return hPed;
}

TMarker *MyRootClass::ButtonFunc3()
{
    return fEventList.at(ip)->GetBaryCenterAsMarker();
}

TLine *MyRootClass::ButtonFunc4()
{
    return fEventList.at(ip)->GetPrincipalAxis();
}

TMarker *MyRootClass::ButtonFunc5()
{
    return fEventList.at(ip)->GetConvertionPoint();
}

TLine *MyRootClass::ButtonFunc6()
{
    return fEventList.at(ip)->GetCovertionAxis();
}

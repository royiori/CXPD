//-----------------------------------------------------
// to do some mathematic tricks for histograms
//
// usage as example:
//
//
//                            Q.LIU
//-----------------------------------------------------

#ifndef MYHISTMATHLIBRARY_H
#define MYHISTMATHLIBRARY_H

#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2F.h>
#include <TH2D.h>

#define SQRTFORM 1
#define SQRTOVERFORM 2

void SmearingGaus(TH1 *hIn, TH1 *hOut, double Res);
void SmearingGaus(TH2 *hIn, TH2 *hOut, double Res);
void SmearingGaus(TH2 *hIn, TH2 *hOut, double Resx, double Resy);
void FastSmearingGaus(TH2 *hIn, TH2 *hOut, double Resx, double Resy);

void SmearingGaus(TH1 *hIn, TH1 *hOut, TF1 *fun);
void SmearingGaus(TH2 *hIn, TH2 *hOut, TF2 *funx, TF2 *funy);
void FastSmearingGaus(TH2 *hIn, TH2 *hOut, TF2 *xFunc, TF2 *yFunc);
void StandardSmearingGaus(TH2 *hIn, TH2 *hOut, TF2 *xFunc, TF2 *yFunc);

//_______Func______________________________________
// Smearing Gaussian 1D
//
void SmearingGaus(TH1 *hIn, TH1 *hOut, double Res)
{
    TF1 *func = new TF1("_func_", Form("%f", Res), -10, 10);
    SmearingGaus(hIn, hOut, func);
    delete func;
}

void SmearingGaus(TH1 *hIn, TH1 *hOut, TF1 *func)
{
    if (hOut == NULL)
        hOut = (TH1 *)hIn->Clone();
    hOut->Reset();

    int nBin = hIn->GetNbinsX();
    double binW = hIn->GetBinWidth(1);
    double *xAxis = new double[nBin];
    double *src = new double[nBin];

    for (int i = 0; i < nBin; i++)
    {
        xAxis[i] = hIn->GetXaxis()->GetBinCenter(i + 1);
        src[i] = hIn->GetBinContent(i + 1);
    }

    for (int i = 0; i < nBin; i++)
    {
        double sigma = func->Eval(xAxis[i]); //sigma may vary for different x-axis loation

        if (src[i] == 0)
            continue;

        if (binW > 6 * sigma)
        {
            hOut->Fill(xAxis[i], src[i]);
            continue;
        }

        for (int j = 0; j < nBin; j++)
        {
            if (fabs(xAxis[j] - xAxis[i]) + binW > 6 * sigma)
                continue;
            double weight = TMath::Erf((xAxis[j] - xAxis[i] + binW / 2) / sigma / sqrt(2)) / 2 - TMath::Erf((xAxis[j] - xAxis[i] - binW / 2) / sigma / sqrt(2)) / 2;
            hOut->Fill(xAxis[j], src[i] * weight);
        }
    }

    delete[] xAxis;
    delete[] src;
}

//_______Func______________________________________
// Smearing Gaussian 2D
//
void SmearingGaus(TH2 *hIn, TH2 *hOut, double Res)
{
    //TF2 *func = new TF2("_func_", Form("%f + 0*x*y", Res), -10, 10, -10, 10);
    FastSmearingGaus(hIn, hOut, Res, Res);
}

void SmearingGaus(TH2 *hIn, TH2 *hOut, double Resx, double Resy)
{
    FastSmearingGaus(hIn, hOut, Resx, Resy);
}

void FastSmearingGaus(TH2 *hIn, TH2 *hOut, double Resx, double Resy)
{
    if (hOut == NULL)
        hOut = (TH2 *)hIn->Clone();
    hOut->Reset();

    // 1. stored to vector
    int nXBin = hIn->GetNbinsX();
    int nYBin = hIn->GetNbinsY();
    double xbinW = hIn->GetXaxis()->GetBinWidth(1);
    double ybinW = hIn->GetYaxis()->GetBinWidth(1);

    double *xAxis = new double[nXBin];
    double *yAxis = new double[nYBin];
    for (int i = 0; i < nXBin; i++)
        xAxis[i] = hIn->GetXaxis()->GetBinCenter(i + 1);
    for (int j = 0; j < nYBin; j++)
        yAxis[j] = hIn->GetYaxis()->GetBinCenter(j + 1);

    std::vector<std::pair<double, double>> Ind;
    std::vector<double> Src;
    for (int i = 0; i < nXBin; i++)
        for (int j = 0; j < nYBin; j++)
        {
            double q = hIn->GetBinContent(i + 1, j + 1);
            if (q == 0)
                continue;

            Ind.push_back(std::make_pair(i, j));
            Src.push_back(q);
        }

    // 2. cal the up/down limit according to smear parameter
    double LimX = int(6 * Resx / xbinW);
    double LimY = int(6 * Resy / ybinW);

    // 3. cal the smear result
    for (unsigned int ind = 0; ind < Ind.size(); ind++)
    {
        
        int i = Ind[ind].first;
        int j = Ind[ind].second;
        int q = Src[ind];

        if (LimX * LimY == 1)
        {
            hOut->SetBinContent(i, j, q);
            continue;
        }

        for (int m = i - LimX; m <= i + LimX; m++)
        {
            for (int n = j - LimY; n < j + LimY; n++)
            {
                if (m < 0 || m >= nXBin)
                    continue;
                if (n < 0 || n >= nYBin)
                    continue;

                double xWeight = TMath::Erf((xAxis[m] - xAxis[i] + xbinW / 2) / Resx / sqrt(2)) / 2 - TMath::Erf((xAxis[m] - xAxis[i] - xbinW / 2) / Resx / sqrt(2)) / 2;
                double yWeight = TMath::Erf((yAxis[n] - yAxis[j] + ybinW / 2) / Resy / sqrt(2)) / 2 - TMath::Erf((yAxis[n] - yAxis[j] - ybinW / 2) / Resy / sqrt(2)) / 2;
                hOut->Fill(xAxis[m], yAxis[n], q * xWeight * yWeight);
            }
        }
    }

    delete[] xAxis;
    delete[] yAxis;
}



void SmearingGaus(TH2 *hIn, TH2 *hOut, TF2 *xFunc, TF2 *yFunc)
{
    StandardSmearingGaus(hIn, hOut, xFunc, yFunc);
}

void StandardSmearingGaus(TH2 *hIn, TH2 *hOut, TF2 *xFunc, TF2 *yFunc)
{
    if (hOut == NULL)
        hOut = (TH2 *)hIn->Clone();
    hOut->Reset();

    int nXBin = hIn->GetNbinsX();
    int nYBin = hIn->GetNbinsY();
    double xbinW = hIn->GetXaxis()->GetBinWidth(1);
    double ybinW = hIn->GetYaxis()->GetBinWidth(1);

    double *xAxis = new double[nXBin];
    double *yAxis = new double[nYBin];

    double *src2 = new double[nXBin * nYBin];
    double **src = new double *[nXBin];
    for (int i = 0; i < nXBin; i++)
        src[i] = src2 + i * nYBin;

    for (int i = 0; i < nXBin; i++)
        xAxis[i] = hIn->GetXaxis()->GetBinCenter(i + 1);
    for (int j = 0; j < nYBin; j++)
        yAxis[j] = hIn->GetYaxis()->GetBinCenter(j + 1);

    for (int i = 0; i < nXBin; i++)
        for (int j = 0; j < nYBin; j++)
            src[i][j] = hIn->GetBinContent(i + 1, j + 1);

    for (int i = 0; i < nXBin; i++)
    {
        for (int j = 0; j < nYBin; j++)
        {
            if (src[i][j] == 0)
                continue;

            double xSigma = xFunc->Eval(xAxis[i], yAxis[j]); //sigma may vary for different x,y-axis loation
            double ySigma = yFunc->Eval(xAxis[i], yAxis[j]); //sigma may vary for different x,y-axis loation

            double xWeight = (xbinW > 6 * xSigma) ? 1. : 0.;
            double yWeight = (ybinW > 6 * ySigma) ? 1. : 0.;

            if (xWeight * yWeight != 0)
            {
                hOut->Fill(xAxis[i], yAxis[j], src[i][j]);
                continue;
            }

            for (int m = 0; m < nXBin; m++)
            {
                for (int n = 0; n < nYBin; n++)
                {
                    xWeight = 1;
                    yWeight = 1;
                    if (fabs(xAxis[m] - xAxis[i]) + xbinW > 6 * xSigma)
                        xWeight = 0.;
                    if (fabs(yAxis[n] - yAxis[j]) + ybinW > 6 * ySigma)
                        yWeight = 0.;

                    if (xWeight * yWeight == 0)
                        continue;

                    xWeight = TMath::Erf((xAxis[m] - xAxis[i] + xbinW / 2) / xSigma / sqrt(2)) / 2 - TMath::Erf((xAxis[m] - xAxis[i] - xbinW / 2) / xSigma / sqrt(2)) / 2;
                    yWeight = TMath::Erf((yAxis[n] - yAxis[j] + ybinW / 2) / ySigma / sqrt(2)) / 2 - TMath::Erf((yAxis[n] - yAxis[j] - ybinW / 2) / ySigma / sqrt(2)) / 2;
                    hOut->Fill(xAxis[m], yAxis[n], src[i][j] * xWeight * yWeight);
                }
            }
        }
    }

    delete[] xAxis;
    delete[] yAxis;
    delete[] src;
}

#endif

#include "MyRootGui.h"

MyRootGui *gMyRootGui;

using namespace std;

int main(int argc, char **argv)
{

    TApplication *theApp;
    theApp = new TApplication("App", &argc, argv);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetStatFont(42);

    gMyRootGui = new MyRootGui(gClient->GetRoot(), 1250, 620); 

    // run ROOT application
    theApp->Run();

    delete theApp;
    return 0;
}
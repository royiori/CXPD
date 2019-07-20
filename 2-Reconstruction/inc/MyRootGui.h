
#ifndef _MyRootGui_h_
#define _MyRootGui_h_

#include "TGNumberEntry.h"
#include "TGComboBox.h"
#include "TGTextEdit.h"
#include "TGSlider.h"
#include "TGFileDialog.h"
#include "TGCanvas.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGMsgBox.h"
#include "TGTab.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TEnv.h"
#include "Riostream.h"
#include "MyRootClass.h"


using namespace std;

//TRootGuiBuilder *g = new TRootGuiBuilder();
class MyRootGui : public TGMainFrame
{
public:
    MyRootGui(const TGWindow *p, int w, int h);
    virtual ~MyRootGui();

    void OpenFile();
    void OpenPedFile();

    void Analysis();
    void DrawPre();
    void DrawNext();
    void DrawSelected();
    void DrawSearchFrameFunc();

    void ShowAllButtonFunc();
    void ShowIPButtonFunc();
    void ShowPolButtonFunc();

    void NSigPedButtonFunc();
    void ShowPedButtonFunc();

    void DataSwitch();

private:
    virtual void CloseWindow();
    
    TGTab *fTab0;
    TGCompositeFrame *fTabPage1;
    TGCompositeFrame *fTabPage2;

    TCanvas *fCanvas1;
    TCanvas *fCanvas2;
    TRootEmbeddedCanvas *fECanvas1;
    TRootEmbeddedCanvas *fECanvas2;
    TGTextButton *fFileButton1;
    TGTextButton *fFileButton2;

    TGTextButton *fAnaButton;
    TGTextButton *fPreButton;
    TGTextButton *fNextButton;
    TGHSlider *fHSlider1;

    TGTextButton *fDrawSFbutton;
    TGTextButton *fShowAllButton;
    TGTextButton *fShowIPButton;
    TGTextButton *fShowPolButton;
    TGTextButton *fNSigPedButton;
    TGTextButton *fShowPedButton;
    TGCheckButton *fDataSwitchButton;
    TGNumberEntry *fNumberEntry31;

    TGLabel *fFileLabel1;
    TGLabel *fFileLabel2;
    TGTextEdit *fOutputText;
    TGTextEdit *fSettingText;

    TString fileDir;
    TString filePath;
    TString pedDir;
    TString pedPath;

    //fTab_page2

    //others
    TEnv *env;
    MyRootClass *gMyRootClass;

    ClassDef(MyRootGui, 0)
};

#endif
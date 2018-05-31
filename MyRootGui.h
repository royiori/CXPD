
#ifndef ROOT_TGDockableFrame
#include "TGDockableFrame.h"
#endif
#ifndef ROOT_TGMdiDecorFrame
#include "TGMdiDecorFrame.h"
#endif
#ifndef ROOT_TG3DLine
#include "TG3DLine.h"
#endif
#ifndef ROOT_TGMdiFrame
#include "TGMdiFrame.h"
#endif
#ifndef ROOT_TGMdiMainFrame
#include "TGMdiMainFrame.h"
#endif
#ifndef ROOT_TGMdiMenu
#include "TGMdiMenu.h"
#endif
#ifndef ROOT_TGListBox
#include "TGListBox.h"
#endif
#ifndef ROOT_TGNumberEntry
#include "TGNumberEntry.h"
#endif
#ifndef ROOT_TGScrollBar
#include "TGScrollBar.h"
#endif
#ifndef ROOT_TGComboBox
#include "TGComboBox.h"
#endif
#ifndef ROOT_TGuiBldHintsEditor
#include "TGuiBldHintsEditor.h"
#endif
#ifndef ROOT_TGuiBldNameFrame
#include "TGuiBldNameFrame.h"
#endif
#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif
#ifndef ROOT_TGMenu
#include "TGMenu.h"
#endif
#ifndef ROOT_TGFileDialog
#include "TGFileDialog.h"
#endif
#ifndef ROOT_TGShutter
#include "TGShutter.h"
#endif
#ifndef ROOT_TGButtonGroup
#include "TGButtonGroup.h"
#endif
#ifndef ROOT_TGCanvas
#include "TGCanvas.h"
#endif
#ifndef ROOT_TGFSContainer
#include "TGFSContainer.h"
#endif
#ifndef ROOT_TGuiBldEditor
#include "TGuiBldEditor.h"
#endif
#ifndef ROOT_TGColorSelect
#include "TGColorSelect.h"
#endif
#ifndef ROOT_TGButton
#include "TGButton.h"
#endif
#ifndef ROOT_TRootContextMenu
#include "TRootContextMenu.h"
#endif
#ifndef ROOT_TGFSComboBox
#include "TGFSComboBox.h"
#endif
#ifndef ROOT_TGLabel
#include "TGLabel.h"
#endif
#ifndef ROOT_TGMsgBox
#include "TGMsgBox.h"
#endif
#ifndef ROOT_TRootGuiBuilder
#include "TRootGuiBuilder.h"
#endif
#ifndef ROOT_TGTab
#include "TGTab.h"
#endif
#ifndef ROOT_TGListView
#include "TGListView.h"
#endif
#ifndef ROOT_TGSplitter
#include "TGSplitter.h"
#endif
#ifndef ROOT_TGStatusBar
#include "TGStatusBar.h"
#endif
#ifndef ROOT_TGListTree
#include "TGListTree.h"
#endif
#ifndef ROOT_TGuiBldGeometryFrame
#include "TGuiBldGeometryFrame.h"
#endif
#ifndef ROOT_TGToolTip
#include "TGToolTip.h"
#endif
#ifndef ROOT_TGToolBar
#include "TGToolBar.h"
#endif
#ifndef ROOT_TRootEmbeddedCanvas
#include "TRootEmbeddedCanvas.h"
#endif
#ifndef ROOT_TCanvas
#include "TCanvas.h"
#endif
#ifndef ROOT_TGuiBldDragManager
#include "TGuiBldDragManager.h"
#endif
#include "TEnv.h"

#include "Riostream.h"

#include "MyRootClass.h"

const char *filetypes[] = {"data files", "*.data",
                           "mdat files", "*.mdat",
                           "ROOT files", "*.root",
                           "All files", "*",
                           0, 0};

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
    void ButtonFunc11();

    void ButtonFunc21();
    void ButtonFunc22();
    void ButtonFunc23();

    void ButtonFunc31();
    void ButtonFunc32();

  private:
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

    TGGroupFrame *fGroupFrame1;
    TGGroupFrame *fGroupFrame2;
    TGGroupFrame *fGroupFrame3;

    TGTextButton *fButton11;
    TGTextButton *fButton21;
    TGTextButton *fButton22;
    TGTextButton *fButton23;
    TGTextButton *fButton31;
    TGTextButton *fButton32;
    TGNumberEntry *fNumberEntry31;

    TGLabel *fFileLabel1;
    TGLabel *fFileLabel2;
    TGLabel *fOutputLabel1;

    TString filePath;
    TString pedPath;

    TEnv *env;
    MyRootClass *gMyRootClass;

    //ClassDef(MyRootGui, 0)
};

//______________________________________________________________________________
//
MyRootGui::MyRootGui(const TGWindow *p, int w, int h) : TGMainFrame(p, w, h)
{
    gStyle->SetOptStat(0);
    env = new TEnv(gSystem->WorkingDirectory()+TString("/.env"));
    env->SaveLevel(kEnvLocal);
    
    filePath = env->GetValue("filePath", filePath);
    pedPath  = env->GetValue("pedPath",  pedPath);

    gMyRootClass = new MyRootClass(filePath, pedPath);

    // embedded canvas
    fECanvas1 = new TRootEmbeddedCanvas(0, this, 500, 500, kSunkenFrame);
    fECanvas1->SetName("fECanvas1");
    Int_t wfECanvas1 = fECanvas1->GetCanvasWindowId();
    TCanvas *c1 = new TCanvas("c1", 10, 10, wfECanvas1);
    fECanvas1->AdoptCanvas(c1);
    this->AddFrame(fECanvas1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fECanvas1->MoveResize(10, 10, 500, 500);
    fCanvas1 = fECanvas1->GetCanvas();

    fFileButton1 = new TGTextButton(this, "Read file", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fFileButton1->Connect("Clicked()", "MyRootGui", this, "OpenFile()");
    fFileButton1->SetTextJustify(36);
    fFileButton1->SetMargins(0, 0, 0, 0);
    fFileButton1->SetWrapLength(-1);
    this->AddFrame(fFileButton1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fFileButton1->MoveResize(10, 520, 98, 23);

    fFileLabel1 = new TGLabel(this, filePath);
    fFileLabel1->SetTextJustify(36);
    fFileLabel1->SetMargins(0, 0, 0, 0);
    fFileLabel1->SetWrapLength(-1);
    fFileLabel1->SetTextJustify(1);
    this->AddFrame(fFileLabel1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fFileLabel1->MoveResize(120, 520, 500, 23);

    fFileButton2 = new TGTextButton(this, "Read ped file", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fFileButton2->Connect("Clicked()", "MyRootGui", this, "OpenPedFile()");
    fFileButton2->SetTextJustify(36);
    fFileButton2->SetMargins(0, 0, 0, 0);
    fFileButton2->SetWrapLength(-1);
    this->AddFrame(fFileButton2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fFileButton2->MoveResize(10, 550, 98, 23);

    fFileLabel2 = new TGLabel(this, pedPath);
    fFileLabel2->SetTextJustify(36);
    fFileLabel2->SetMargins(0, 0, 0, 0);
    fFileLabel2->SetWrapLength(-1);
    fFileLabel2->SetTextJustify(1);
    this->AddFrame(fFileLabel2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fFileLabel2->MoveResize(120, 550, 500, 23);

    //---------------------
    // Button lists:
    double x0, y0, x, y;

    x0 = 520;
    y0 = 20;

    fAnaButton = new TGTextButton(this, "Analysis", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fAnaButton->Connect("Clicked()", "MyRootGui", this, "Analysis()");
    fAnaButton->SetTextJustify(36);
    fAnaButton->SetMargins(0, 0, 0, 0);
    fAnaButton->SetWrapLength(-1);
    this->AddFrame(fAnaButton, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fAnaButton->MoveResize(x0, y0, 220, 23);
    fAnaButton->SetEnabled(kTRUE);

    //---------------------
    // "fGroupFrame1" group frame
    x0 = 520;
    y0 = 60;

    fGroupFrame1 = new TGGroupFrame(this, "Display single event:");
    fGroupFrame1->SetLayoutBroken(kTRUE);
    fGroupFrame1->SetLayoutManager(new TGVerticalLayout(fGroupFrame1));
    this->AddFrame(fGroupFrame1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fGroupFrame1->MoveResize(x0, y0, 220, 105);

    fHSlider1 = new TGHSlider(this,134,kSlider1 | kScaleBoth,-1,kHorizontalFrame);
    fHSlider1->Connect("Released()", "MyRootGui", this, "DrawSelected()");
    fHSlider1->SetName("fHSlider1");
    fHSlider1->SetRange(0,40);
    fHSlider1->SetPosition(20);
    this->AddFrame(fHSlider1, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fHSlider1->MoveResize(x0 + 10, y0 + 20, 200, 24);
    fHSlider1->SetEnabled(kFALSE);

    fPreButton = new TGTextButton(this, "<", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fPreButton->Connect("Clicked()", "MyRootGui", this, "DrawPre()");
    fPreButton->SetTextJustify(36);
    fPreButton->SetMargins(0, 0, 0, 0);
    fPreButton->SetWrapLength(-1);
    fGroupFrame1->AddFrame(fPreButton, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fPreButton->MoveResize(x0 + 10, y0 + 45, 100, 23);
    fPreButton->SetEnabled(kFALSE);

    fNextButton = new TGTextButton(this, ">", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fNextButton->Connect("Clicked()", "MyRootGui", this, "DrawNext()");
    fNextButton->SetTextJustify(36);
    fNextButton->SetMargins(0, 0, 0, 0);
    fNextButton->SetWrapLength(-1);
    fGroupFrame1->AddFrame(fNextButton, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fNextButton->MoveResize(x0 + 110, y0 + 45, 100, 23);
    fNextButton->SetEnabled(kFALSE);

    fButton11 = new TGTextButton(this, "Draw search frame", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton11->Connect("Clicked()", "MyRootGui", this, "ButtonFunc11()");
    fButton11->SetTextJustify(36);
    fButton11->SetMargins(0, 0, 0, 0);
    fButton11->SetWrapLength(-1);
    this->AddFrame(fButton11, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton11->MoveResize(x0 + 10, y0 + 70, 200, 23);
    fButton11->SetEnabled(kFALSE);

    //---------------------
    // "fGroupFrame2" group frame
    x0 = 520;
    y0 = 175;

    fGroupFrame2 = new TGGroupFrame(this, "Event summary: ");
    fGroupFrame2->SetLayoutBroken(kTRUE);
    fGroupFrame2->SetLayoutManager(new TGVerticalLayout(fGroupFrame2));
    this->AddFrame(fGroupFrame2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fGroupFrame2->MoveResize(x0, y0, 220, 105);

    fButton21 = new TGTextButton(this, "Show all events", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton21->Connect("Clicked()", "MyRootGui", this, "ButtonFunc21()");
    fButton21->SetTextJustify(36);
    fButton21->SetMargins(0, 0, 0, 0);
    fButton21->SetWrapLength(-1);
    this->AddFrame(fButton21, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton21->MoveResize(x0 + 10, y0 + 20, 200, 23);
    fButton21->SetEnabled(kFALSE);

    fButton22 = new TGTextButton(this, "Show center/IP points", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton22->Connect("Clicked()", "MyRootGui", this, "ButtonFunc22()");
    fButton22->SetTextJustify(36);
    fButton22->SetMargins(0, 0, 0, 0);
    fButton22->SetWrapLength(-1);
    this->AddFrame(fButton22, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton22->MoveResize(x0 + 10, y0 + 45, 200, 23);
    fButton22->SetEnabled(kFALSE);

    fButton23 = new TGTextButton(this, "Show polarization results", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton23->Connect("Clicked()", "MyRootGui", this, "ButtonFunc23()");
    fButton23->SetTextJustify(36);
    fButton23->SetMargins(0, 0, 0, 0);
    fButton23->SetWrapLength(-1);
    this->AddFrame(fButton23, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton23->MoveResize(x0 + 10, y0 + 70, 200, 23);
    fButton23->SetEnabled(kFALSE);

    //---------------------
    // "fGroupFrame3" group frame
    x0 = 520;
    y0 = 290;

    fGroupFrame3 = new TGGroupFrame(this, "Pedestal: ");
    fGroupFrame3->SetLayoutBroken(kTRUE);
    fGroupFrame3->SetLayoutManager(new TGVerticalLayout(fGroupFrame2));
    this->AddFrame(fGroupFrame2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fGroupFrame3->MoveResize(x0, y0, 220, 80);

    fButton31 = new TGTextButton(this, "Analysis Pedestal", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton31->Connect("Clicked()", "MyRootGui", this, "ButtonFunc31()");
    fButton31->SetTextJustify(36);
    fButton31->SetMargins(0, 0, 0, 0);
    fButton31->SetWrapLength(-1);
    this->AddFrame(fButton31, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton31->MoveResize(x0 + 10, y0 + 20, 160, 23);
    fButton31->SetEnabled(kTRUE);

    fNumberEntry31 = new TGNumberEntry(this, (Double_t)3, 6, -1, (TGNumberFormat::EStyle)5);
    fNumberEntry31->SetName("fNumberEntry31");
    this->AddFrame(fNumberEntry31, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fNumberEntry31->MoveResize(x0 + 170, y0 + 20, 40, 23);

    fButton32 = new TGTextButton(this, "Show pedstal m/s", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton32->Connect("Clicked()", "MyRootGui", this, "ButtonFunc32()");
    fButton32->SetTextJustify(36);
    fButton32->SetMargins(0, 0, 0, 0);
    fButton32->SetWrapLength(-1);
    this->AddFrame(fButton32, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton32->MoveResize(x0 + 10, y0 + 45, 200, 23);
    fButton32->SetEnabled(kFALSE);

    //---------------------
    // "fGroupFrame4" group frame

    x0 = 520;
    y0 = 380;

    fOutputLabel1 = new TGLabel(this, "");
    fOutputLabel1->SetTextJustify(36);
    fOutputLabel1->SetMargins(0, 0, 0, 0);
    fOutputLabel1->SetWrapLength(-1);
    fOutputLabel1->SetTextJustify(1);
    this->AddFrame(fOutputLabel1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fOutputLabel1->MoveResize(x0, y0, 220, 200);

    //---------------------
    // embedded canvas
    fECanvas2 = new TRootEmbeddedCanvas(0, this, 500, 500, kSunkenFrame);
    fECanvas2->SetName("fECanvas2");
    Int_t wfECanvas2 = fECanvas2->GetCanvasWindowId();
    TCanvas *c2 = new TCanvas("c2", 10, 10, wfECanvas2);
    fECanvas2->AdoptCanvas(c2);
    this->AddFrame(fECanvas2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fECanvas2->MoveResize(750, 10, 500, 500);
    fCanvas2 = fECanvas2->GetCanvas();

    this->MapSubwindows();
    this->MapWindow();
}

MyRootGui::~MyRootGui()
{
    // Destructor.
}

//______________________________________________________________________________
//
void MyRootGui::OpenFile()
{
    static TString dir(".");
    TGFileInfo fi;
    fi.fFileTypes = filetypes;
    fi.fIniDir = StrDup(dir);
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if(fi.fFilename == NULL) return;
    filePath = fi.fFilename;
    fFileLabel1->SetText(filePath);
    env->SetValue("filePath", filePath);
    env->Save();

    fAnaButton->SetEnabled(kTRUE);
    fHSlider1->SetEnabled(kFALSE);
    fPreButton->SetEnabled(kFALSE);
    fNextButton->SetEnabled(kFALSE);
    fButton11->SetEnabled(kFALSE);

    fButton21->SetEnabled(kFALSE);
    fButton22->SetEnabled(kFALSE);
    fButton23->SetEnabled(kFALSE);
}

void MyRootGui::OpenPedFile()
{
    static TString dir(".");
    TGFileInfo fi;
    fi.fFileTypes = filetypes;
    fi.fIniDir = StrDup(dir);
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if(fi.fFilename == NULL) return;    
    pedPath = fi.fFilename;
    fFileLabel2->SetText(pedPath);
    env->SetValue("pedPath",  pedPath);
    env->Save();

    fButton31->SetEnabled(kTRUE);
    fButton32->SetEnabled(kFALSE);
}

//______________________________________________________________________________
//
void MyRootGui::Analysis()
{
    fCanvas1->cd();
    gMyRootClass->Analysis(filePath);
    fAnaButton->SetEnabled(kFALSE);
    fHSlider1->SetEnabled(kTRUE); 
    fPreButton->SetEnabled(kTRUE);
    fNextButton->SetEnabled(kTRUE);
    fButton11->SetEnabled(kTRUE);

    fButton21->SetEnabled(kTRUE);
    fButton22->SetEnabled(kTRUE);
    fButton23->SetEnabled(kTRUE);

    fHSlider1->SetRange(0, gMyRootClass->GetNEvent());
    fHSlider1->SetPosition(gMyRootClass->GetIp());

    ButtonFunc21();
}

void MyRootGui::DrawPre()
{
    fCanvas1->cd();
    gMyRootClass->DrawPre("box");
    fCanvas1->Modified();
    fCanvas1->Update();
    fOutputLabel1->SetText(gMyRootClass->GetInfo()->Data());
    fHSlider1->SetPosition(gMyRootClass->GetIp());

    fCanvas2->cd();
    gMyRootClass->DrawRaw();
    fCanvas2->Modified();
    fCanvas2->Update();
}

void MyRootGui::DrawNext()
{
    fCanvas1->cd();
    gMyRootClass->DrawNext("box");
    fCanvas1->Modified();
    fCanvas1->Update();
    fOutputLabel1->SetText(gMyRootClass->GetInfo()->Data());
    fHSlider1->SetPosition(gMyRootClass->GetIp());

    fCanvas2->cd();
    gMyRootClass->DrawRaw();
    fCanvas2->Modified();
    fCanvas2->Update();
}

void MyRootGui::DrawSelected()
{
    fCanvas1->cd();
    gMyRootClass->DrawSelected(fHSlider1->GetPosition(), "box");
    fCanvas1->Modified();
    fCanvas1->Update();
    fOutputLabel1->SetText(gMyRootClass->GetInfo()->Data());

    fCanvas2->cd();
    gMyRootClass->DrawRaw();
    fCanvas2->Modified();
    fCanvas2->Update();
}

void MyRootGui::ButtonFunc11()
{
    fCanvas1->cd();
    gMyRootClass->ButtonFunc11();
    fCanvas1->Modified();
    fCanvas1->Update();
}

//______________________________________________________________________________
//
void MyRootGui::ButtonFunc21()
{
    fCanvas1->cd();
    gMyRootClass->ButtonFunc21(1);
    fCanvas1->Modified();
    fCanvas1->Update();

    fCanvas2->cd();
    gMyRootClass->ButtonFunc21(2, "box");
    fCanvas2->Modified();
    fCanvas2->Update();    
}

void MyRootGui::ButtonFunc22()
{
    fCanvas1->cd();
    gMyRootClass->ButtonFunc22(1);
    fCanvas1->Modified();
    fCanvas1->Update();

    fCanvas2->cd();
    gMyRootClass->ButtonFunc22(2);
    fCanvas2->Modified();
    fCanvas2->Update();        
}

void MyRootGui::ButtonFunc23()
{
    fCanvas1->cd();
    gMyRootClass->ButtonFunc23(1);
    fCanvas1->Modified();
    fCanvas1->Update();

    fCanvas2->cd();
    gMyRootClass->ButtonFunc23(2);
    fCanvas2->Modified();
    fCanvas2->Update();    
}

//______________________________________________________________________________

//______________________________________________________________________________
//
void MyRootGui::ButtonFunc31()
{
    gMyRootClass->ButtonFunc31(fNumberEntry31->GetNumber());

    fCanvas1->cd();
    gMyRootClass->AnalysisPed(pedPath);
    fCanvas1->Modified();
    fCanvas1->Update();

    fCanvas2->cd();
    gMyRootClass->ButtonFunc33();
    fCanvas2->Modified();
    fCanvas2->Update();
        
    fButton31->SetEnabled(kFALSE);
    fButton32->SetEnabled(kTRUE);
}

void MyRootGui::ButtonFunc32()
{
    gMyRootClass->ButtonFunc31(fNumberEntry31->GetNumber());

    fCanvas1->cd();
    gMyRootClass->ButtonFunc32();
    fCanvas1->Modified();
    fCanvas1->Update();

    fCanvas2->cd();
    gMyRootClass->ButtonFunc33();
    fCanvas2->Modified();
    fCanvas2->Update();
}

//______________________________________________________________________________

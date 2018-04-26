
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

#include "Riostream.h"

#include "MyRootClass.h"

const char *filetypes[] = {"mdat files", "*.mdat",
                           "ROOT files", "*.root",
                           "All files", "*",
                           0, 0};

class MyRootGui : public TGMainFrame
{
  public:
    MyRootGui(const TGWindow *p, int w, int h);
    virtual ~MyRootGui();

    void OpenFile();
    void Analysis();
    void DrawPre();
    void DrawNext();

    void ButtonFunc1();
    void ButtonFunc2();
    void ButtonFunc3();
    void ButtonFunc4();
    void ButtonFunc5();
    void ButtonFunc6();

  private:
    TCanvas *fCanvas1;
    TRootEmbeddedCanvas *fECanvas1;
    TGTextButton *fFileButton1;
    TGTextButton *fAnaButton;
    TGTextButton *fPreButton;
    TGTextButton *fNextButton;

    TGTextButton *fButton1;
    TGTextButton *fButton2;
    TGTextButton *fButton3;
    TGTextButton *fButton4;
    TGTextButton *fButton5;
    TGTextButton *fButton6;
    TGNumberEntry *fNumberEntry2;

    TGLabel *fFileLabel1;
    TGLabel *fOutputLabel1;

    TString filePath;

    MyRootClass *gMyRootClass;

    //ClassDef(MyRootGui, 0)
};

//______________________________________________________________________________
//
MyRootGui::MyRootGui(const TGWindow *p, int w, int h) : TGMainFrame(p, w, h)
{
    filePath = "/Users/liuqian/Work/work/CXPD/pixy/dme_simul/cxpd.mdat";

    gMyRootClass = new MyRootClass();

    // embedded canvas
    fECanvas1 = new TRootEmbeddedCanvas(0, this, 600, 600, kSunkenFrame);
    fECanvas1->SetName("fECanvas1");
    Int_t wfECanvas1 = fECanvas1->GetCanvasWindowId();
    TCanvas *c123 = new TCanvas("c123", 10, 10, wfECanvas1);
    fECanvas1->AdoptCanvas(c123);
    this->AddFrame(fECanvas1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fECanvas1->MoveResize(10, 10, 600, 600);
    fCanvas1 = fECanvas1->GetCanvas();

    fFileButton1 = new TGTextButton(this, "Read file", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fFileButton1->Connect("Clicked()", "MyRootGui", this, "OpenFile()");
    fFileButton1->SetTextJustify(36);
    fFileButton1->SetMargins(0, 0, 0, 0);
    fFileButton1->SetWrapLength(-1);
    fFileButton1->Resize(98, 23);
    this->AddFrame(fFileButton1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fFileButton1->MoveResize(10, 620, 98, 23);

    fFileLabel1 = new TGLabel(this, filePath);
    fFileLabel1->SetTextJustify(36);
    fFileLabel1->SetMargins(0, 0, 0, 0);
    fFileLabel1->SetWrapLength(-1);
    fFileLabel1->SetTextJustify(1);
    this->AddFrame(fFileLabel1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fFileLabel1->MoveResize(120, 622, 500, 23);

    fAnaButton = new TGTextButton(this, "Analysis", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fAnaButton->Connect("Clicked()", "MyRootGui", this, "Analysis()");
    fAnaButton->SetTextJustify(36);
    fAnaButton->SetMargins(0, 0, 0, 0);
    fAnaButton->SetWrapLength(-1);
    fAnaButton->Resize(98, 23);
    this->AddFrame(fAnaButton, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fAnaButton->MoveResize(620, 20, 220, 23);

    fPreButton = new TGTextButton(this, "<", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fPreButton->Connect("Clicked()", "MyRootGui", this, "DrawPre()");
    fPreButton->SetTextJustify(36);
    fPreButton->SetMargins(0, 0, 0, 0);
    fPreButton->SetWrapLength(-1);
    fPreButton->Resize(98, 23);
    this->AddFrame(fPreButton, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fPreButton->MoveResize(620, 60, 110, 23);

    fNextButton = new TGTextButton(this, ">", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fNextButton->Connect("Clicked()", "MyRootGui", this, "DrawNext()");
    fNextButton->SetTextJustify(36);
    fNextButton->SetMargins(0, 0, 0, 0);
    fNextButton->SetWrapLength(-1);
    fNextButton->Resize(98, 23);
    this->AddFrame(fNextButton, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fNextButton->MoveResize(730, 60, 110, 23);

    fButton1 = new TGTextButton(this, "Calculate Pedestal", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton1->Connect("Clicked()", "MyRootGui", this, "ButtonFunc1()");
    fButton1->SetTextJustify(36);
    fButton1->SetMargins(0, 0, 0, 0);
    fButton1->SetWrapLength(-1);
    fButton1->Resize(98, 23);
    this->AddFrame(fButton1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton1->MoveResize(620, 100, 220, 23);

    fButton2 = new TGTextButton(this, "Set Pedestal:", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton2->Connect("Clicked()", "MyRootGui", this, "ButtonFunc2()");
    fButton2->SetTextJustify(36);
    fButton2->SetMargins(0, 0, 0, 0);
    fButton2->SetWrapLength(-1);
    fButton2->Resize(98, 23);
    this->AddFrame(fButton2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton2->MoveResize(620, 140, 90, 23);

    fNumberEntry2 = new TGNumberEntry(this, (Double_t)0, 6, -1, (TGNumberFormat::EStyle)5);
    fNumberEntry2->SetName("fNumberEntry2");
    this->AddFrame(fNumberEntry2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fNumberEntry2->MoveResize(720, 140, 120, 23);

    fButton3 = new TGTextButton(this, "Find Barycenter", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton3->Connect("Clicked()", "MyRootGui", this, "ButtonFunc3()");
    fButton3->SetTextJustify(36);
    fButton3->SetMargins(0, 0, 0, 0);
    fButton3->SetWrapLength(-1);
    fButton3->Resize(98, 23);
    this->AddFrame(fButton3, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton3->MoveResize(620, 180, 220, 23);

    fButton4 = new TGTextButton(this, "Recon principal axis", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton4->Connect("Clicked()", "MyRootGui", this, "ButtonFunc4()");
    fButton4->SetTextJustify(36);
    fButton4->SetMargins(0, 0, 0, 0);
    fButton4->SetWrapLength(-1);
    fButton4->Resize(98, 23);
    this->AddFrame(fButton4, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton4->MoveResize(620, 205, 220, 23);

    fButton5 = new TGTextButton(this, "Find Convertion Point", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton5->Connect("Clicked()", "MyRootGui", this, "ButtonFunc5()");
    fButton5->SetTextJustify(36);
    fButton5->SetMargins(0, 0, 0, 0);
    fButton5->SetWrapLength(-1);
    fButton5->Resize(98, 23);
    this->AddFrame(fButton5, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton5->MoveResize(620, 230, 220, 23);

    fButton6 = new TGTextButton(this, "Recon Incident axis",  -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
    fButton6->Connect("Clicked()", "MyRootGui", this, "ButtonFunc6()");
    fButton6->SetTextJustify(36);
    fButton6->SetMargins(0, 0, 0, 0);
    fButton6->SetWrapLength(-1);
    fButton6->Resize(98, 23);
    this->AddFrame(fButton6, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fButton6->MoveResize(620, 255, 220, 23);

    fOutputLabel1 = new TGLabel(this, "");
    fOutputLabel1->SetTextJustify(36);
    fOutputLabel1->SetMargins(0, 0, 0, 0);
    fOutputLabel1->SetWrapLength(-1);
    fOutputLabel1->SetTextJustify(1);
    this->AddFrame(fOutputLabel1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fOutputLabel1->MoveResize(630, 290, 220, 350);

    this->MapSubwindows();
    this->MapWindow();
}

MyRootGui::~MyRootGui()
{
    // Destructor.
}

void MyRootGui::OpenFile()
{
    static TString dir(".");
    TGFileInfo fi;
    fi.fFileTypes = filetypes;
    fi.fIniDir = StrDup(dir);
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    filePath = fi.fFilename;
    fFileLabel1->SetText(filePath);
}

void MyRootGui::Analysis()
{
    gMyRootClass->Analysis(filePath);
    DrawNext();
}


void MyRootGui::DrawPre()
{
    TH2F *f = gMyRootClass->DrawPre();
    if(f==NULL) return;
    f->Draw("colz");
    fCanvas1->Modified();
    fCanvas1->Update();
    fOutputLabel1->SetText(gMyRootClass->GetInfo()->Data());
    ButtonFunc3();
    ButtonFunc4();
}

void MyRootGui::DrawNext()
{
    TH2F *f = gMyRootClass->DrawNext();
    if(f==NULL) return;
    f->Draw("colz");
    fCanvas1->Modified();
    fCanvas1->Update();
    fOutputLabel1->SetText(gMyRootClass->GetInfo()->Data());
    ButtonFunc3();
    ButtonFunc4();
    ButtonFunc5();
    ButtonFunc6();
}

void MyRootGui::ButtonFunc1()
{
    TH1F *f = gMyRootClass->ButtonFunc1();
    if(f==NULL) return;
    f->Draw();
    fCanvas1->Modified();
    fCanvas1->Update();
}

void MyRootGui::ButtonFunc2()
{
    gMyRootClass->ButtonFunc2(fNumberEntry2->GetNumber());
}

void MyRootGui::ButtonFunc3()
{
    TMarker *m = gMyRootClass->ButtonFunc3();
    if(m==NULL) return;
    m->Draw("same");
    fCanvas1->Modified();
    fCanvas1->Update();
}

void MyRootGui::ButtonFunc4()
{
    TLine *l = gMyRootClass->ButtonFunc4();
    if(l==NULL) return;
    l->Draw("same");
    fCanvas1->Modified();
    fCanvas1->Update();
}

void MyRootGui::ButtonFunc5()
{
    TMarker *m = gMyRootClass->ButtonFunc5();
    if(m==NULL) return;
    m->Draw("same");
    fCanvas1->Modified();
    fCanvas1->Update();
}

void MyRootGui::ButtonFunc6()
{
    TLine *l = gMyRootClass->ButtonFunc6();
    if(l==NULL) return;
    l->Draw("same");
    fCanvas1->Modified();
    fCanvas1->Update();
}
//______________________________________________________________________________

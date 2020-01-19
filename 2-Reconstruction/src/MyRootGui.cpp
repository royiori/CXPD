
#include "MyRootGui.h"

using namespace std;

const char *filetypes[] = {"data files", "*.data",
                           "mdat files", "*.mdat",
                           "ROOT files", "*.root",
                           "All files", "*",
                           0, 0};

//______________________________________________________________________________
//
MyRootGui::MyRootGui(const TGWindow *p, int w, int h) : TGMainFrame(p, w, h)
{
    gStyle->SetOptStat(0);
    env = new TEnv(gSystem->WorkingDirectory() + TString("/.env"));
    env->SaveLevel(kEnvLocal);

    fileDir = env->GetValue("fileDir", "");
    pedDir = env->GetValue("pedDir", "");
    filePath = env->GetValue("filePath", "");
    pedPath = env->GetValue("pedPath", "");
    int _dataswitch = env->GetValue("DataSwitch", 0);

    gMyRootClass = new MyRootClass(filePath, pedPath);
    gMyRootClass->UpdateEnvParameters(env);

    // tab
    fTab0 = new TGTab(this, 10, 10);

    // tab0_page1
    fTabPage1 = fTab0->AddTab("Analysis");
    fTabPage1->SetLayoutManager(new TGVerticalLayout(fTabPage1));

    //fillings in tab0_page1
    TGVerticalFrame *fTabPage1Frame = new TGVerticalFrame(fTabPage1, 10, 10, kVerticalFrame);
    {

        //fillings in fTabPage1Frame: ResultFrame + ReadButtonFrame + ReadPadFrame
        //-->ResultFrame
        TGHorizontalFrame *ResultFrame = new TGHorizontalFrame(fTabPage1Frame, 10, 10, kHorizontalFrame);
        {
            //1. embedded canvas
            fECanvas1 = new TRootEmbeddedCanvas(0, ResultFrame, 500, 500, kSunkenFrame);
            fECanvas1->SetName("fECanvas1");
            Int_t wfECanvas1 = fECanvas1->GetCanvasWindowId();
            TCanvas *c1 = new TCanvas("c1", 10, 10, wfECanvas1);
            fECanvas1->AdoptCanvas(c1);
            ResultFrame->AddFrame(fECanvas1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY, 2, 2, 2, 2));
            fCanvas1 = fECanvas1->GetCanvas();

            //2. result button lists
            TGVerticalFrame *ResultButtonFrame = new TGVerticalFrame(ResultFrame, 220, 100, kVerticalFrame);
            {
                fAnaButton = new TGTextButton(ResultButtonFrame, "Analysis", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                fAnaButton->Connect("Clicked()", "MyRootGui", this, "Analysis()");
                fAnaButton->SetTextJustify(36);
                fAnaButton->SetMargins(0, 0, 0, 0);
                fAnaButton->SetWrapLength(-1);
                ResultButtonFrame->AddFrame(fAnaButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                fAnaButton->SetEnabled(kTRUE);

                TGGroupFrame *fGroupFrame1 = new TGGroupFrame(ResultButtonFrame, "Display single event:");
                //fGroupFrame1->SetLayoutBroken(kTRUE);
                fGroupFrame1->SetLayoutManager(new TGVerticalLayout(fGroupFrame1));
                {
                    fHSlider1 = new TGHSlider(fGroupFrame1, 200, kSlider1 | kScaleBoth, -1, kHorizontalFrame);
                    fHSlider1->Connect("Released()", "MyRootGui", this, "DrawSelected()");
                    fHSlider1->SetName("fHSlider1");
                    fHSlider1->SetRange(0, 40);
                    fHSlider1->SetPosition(20);
                    fGroupFrame1->AddFrame(fHSlider1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                    fHSlider1->SetEnabled(kFALSE);

                    TGHorizontalFrame *fSliderFrame = new TGHorizontalFrame(fGroupFrame1, 220, 100, kHorizontalFrame);
                    {
                        fPreButton = new TGTextButton(fSliderFrame, "<", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                        fPreButton->Connect("Clicked()", "MyRootGui", this, "DrawPre()");
                        fPreButton->SetTextJustify(36);
                        fPreButton->SetMargins(0, 0, 0, 0);
                        fPreButton->SetWrapLength(-1);
                        fSliderFrame->AddFrame(fPreButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                        fPreButton->SetEnabled(kFALSE);

                        fNextButton = new TGTextButton(fSliderFrame, ">", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                        fNextButton->Connect("Clicked()", "MyRootGui", this, "DrawNext()");
                        fNextButton->SetTextJustify(36);
                        fNextButton->SetMargins(0, 0, 0, 0);
                        fNextButton->SetWrapLength(-1);
                        fSliderFrame->AddFrame(fNextButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));
                        fNextButton->SetEnabled(kFALSE);
                    }
                    fGroupFrame1->AddFrame(fSliderFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));

                    fDrawSFbutton = new TGTextButton(fGroupFrame1, "Draw search result", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                    fDrawSFbutton->Connect("Clicked()", "MyRootGui", this, "DrawSearchFrameFunc()");
                    fDrawSFbutton->SetTextJustify(36);
                    fDrawSFbutton->SetMargins(0, 0, 0, 0);
                    fDrawSFbutton->SetWrapLength(-1);
                    fGroupFrame1->AddFrame(fDrawSFbutton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                    fDrawSFbutton->SetEnabled(kFALSE);
                }
                ResultButtonFrame->AddFrame(fGroupFrame1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));

                TGGroupFrame *fGroupFrame2 = new TGGroupFrame(ResultButtonFrame, "Event summary: ");
                fGroupFrame2->SetLayoutManager(new TGVerticalLayout(fGroupFrame2));
                {
                    fShowAllButton = new TGTextButton(fGroupFrame2, "Show all events", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                    fShowAllButton->Connect("Clicked()", "MyRootGui", this, "ShowAllButtonFunc()");
                    fShowAllButton->SetTextJustify(36);
                    fShowAllButton->SetMargins(0, 0, 0, 0);
                    fShowAllButton->SetWrapLength(-1);
                    fGroupFrame2->AddFrame(fShowAllButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                    fShowAllButton->SetEnabled(kFALSE);

                    fShowIPButton = new TGTextButton(fGroupFrame2, "Show center/IP points", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                    fShowIPButton->Connect("Clicked()", "MyRootGui", this, "ShowIPButtonFunc()");
                    fShowIPButton->SetTextJustify(36);
                    fShowIPButton->SetMargins(0, 0, 0, 0);
                    fShowIPButton->SetWrapLength(-1);
                    fGroupFrame2->AddFrame(fShowIPButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                    fShowIPButton->SetEnabled(kFALSE);

                    fShowPolButton = new TGTextButton(fGroupFrame2, "Show polarization results", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                    fShowPolButton->Connect("Clicked()", "MyRootGui", this, "ShowPolButtonFunc()");
                    fShowPolButton->SetTextJustify(36);
                    fShowPolButton->SetMargins(0, 0, 0, 0);
                    fShowPolButton->SetWrapLength(-1);
                    fGroupFrame2->AddFrame(fShowPolButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                    fShowPolButton->SetEnabled(kFALSE);
                }
                ResultButtonFrame->AddFrame(fGroupFrame2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));

                TGGroupFrame *fGroupFrame3 = new TGGroupFrame(ResultButtonFrame, "Pedestal: ");
                fGroupFrame3->SetLayoutManager(new TGVerticalLayout(fGroupFrame3));
                {
                    fDataSwitchButton = new TGCheckButton(fGroupFrame3, "Use Pedstal?");
                    fDataSwitchButton->SetTextJustify(36);
                    fDataSwitchButton->SetMargins(0, 0, 0, 0);
                    fDataSwitchButton->SetWrapLength(-1);
                    fDataSwitchButton->Connect("Clicked()", "MyRootGui", this, "DataSwitch()");
                    fGroupFrame3->AddFrame(fDataSwitchButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));

                    fShowPedButton = new TGTextButton(fGroupFrame3, "Show pedstal", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                    fShowPedButton->Connect("Clicked()", "MyRootGui", this, "ShowPedButtonFunc()");
                    fShowPedButton->SetTextJustify(36);
                    fShowPedButton->SetMargins(0, 0, 0, 0);
                    fShowPedButton->SetWrapLength(-1);
                    fGroupFrame3->AddFrame(fShowPedButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                    fShowPedButton->SetEnabled(kFALSE);

                    TGHorizontalFrame *fPedFrame2 = new TGHorizontalFrame(fGroupFrame3, 10, 10, kHorizontalFrame);
                    {
                        fNSigPedButton = new TGTextButton(fPedFrame2, "Analysis Pedestal", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                        fNSigPedButton->Connect("Clicked()", "MyRootGui", this, "NSigPedButtonFunc()");
                        fNSigPedButton->SetTextJustify(36);
                        fNSigPedButton->SetMargins(0, 0, 0, 0);
                        fNSigPedButton->SetWrapLength(-1);
                        fPedFrame2->AddFrame(fNSigPedButton, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                        fNSigPedButton->SetEnabled(kFALSE);

                        fNumberEntry31 = new TGNumberEntry(fPedFrame2, (Double_t)3, 6, -1, (TGNumberFormat::EStyle)5);
                        fNumberEntry31->SetName("fNumberEntry31");
                        fPedFrame2->AddFrame(fNumberEntry31, new TGLayoutHints(kLHintsLeft | kLHintsTop | kFixedWidth, 2, 2, 2, 2));
                        fNumberEntry31->SetState(kFALSE);
                    }
                    fGroupFrame3->AddFrame(fPedFrame2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
                }
                ResultButtonFrame->AddFrame(fGroupFrame3, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));

                fOutputText = new TGTextEdit(ResultButtonFrame, 220, 140);
                fOutputText->LoadBuffer("------> Ready <-------");
                ResultButtonFrame->AddFrame(fOutputText, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
            }
            ResultFrame->AddFrame(ResultButtonFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

            //3. embedded canvas
            fECanvas2 = new TRootEmbeddedCanvas(0, ResultFrame, 500, 500, kSunkenFrame);
            fECanvas2->SetName("fECanvas2");
            Int_t wfECanvas2 = fECanvas2->GetCanvasWindowId();
            TCanvas *c2 = new TCanvas("c2", 10, 10, wfECanvas2);
            fECanvas2->AdoptCanvas(c2);
            ResultFrame->AddFrame(fECanvas2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY, 2, 2, 2, 2));
            fECanvas2->MoveResize(750, 10, 500, 500);
            fCanvas2 = fECanvas2->GetCanvas();
        }
        fTabPage1Frame->AddFrame(ResultFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

        //-->ReadStatFrame
        TGHorizontalFrame *ReadStatFrame = new TGHorizontalFrame(fTabPage1Frame, 10, 10, kHorizontalFrame);
        {
            //-->ReadButtonFrame
            TGVerticalFrame *ReadButtonFrame = new TGVerticalFrame(ReadStatFrame, 98, 48, kFixedWidth);
            {
                fFileButton1 = new TGTextButton(ReadButtonFrame, "Read file", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                fFileButton1->Connect("Clicked()", "MyRootGui", this, "OpenFile()");
                fFileButton1->SetTextJustify(36);
                fFileButton1->SetMargins(0, 0, 0, 0);
                fFileButton1->SetWrapLength(-1);
                ReadButtonFrame->AddFrame(fFileButton1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));

                fFileButton2 = new TGTextButton(ReadButtonFrame, "Read ped file", -1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
                fFileButton2->Connect("Clicked()", "MyRootGui", this, "OpenPedFile()");
                fFileButton2->SetTextJustify(36);
                fFileButton2->SetMargins(0, 0, 0, 0);
                fFileButton2->SetWrapLength(-1);
                ReadButtonFrame->AddFrame(fFileButton2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
            }
            ReadStatFrame->AddFrame(ReadButtonFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));

            //-->ReadLabelFrame
            TGVerticalFrame *ReadLabelFrame = new TGVerticalFrame(ReadStatFrame, 600, 48, kFixedWidth);
            {
                fFileLabel1 = new TGLabel(ReadLabelFrame, filePath);
                fFileLabel1->SetTextJustify(36);
                fFileLabel1->SetMargins(0, 0, 0, 0);
                fFileLabel1->SetWrapLength(-1);
                fFileLabel1->SetTextJustify(1);
                ReadLabelFrame->AddFrame(fFileLabel1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 7, 2));

                fFileLabel2 = new TGLabel(ReadLabelFrame, pedPath);
                fFileLabel2->SetTextJustify(36);
                fFileLabel2->SetMargins(0, 0, 0, 0);
                fFileLabel2->SetWrapLength(-1);
                fFileLabel2->SetTextJustify(1);
                ReadLabelFrame->AddFrame(fFileLabel2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 7, 2));
            }
            ReadStatFrame->AddFrame(ReadLabelFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
        }
        fTabPage1Frame->AddFrame(ReadStatFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
    }
    fTabPage1->AddFrame(fTabPage1Frame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    // tab0_page2
    fTabPage2 = fTab0->AddTab("Settings");
    fTabPage2->SetLayoutManager(new TGVerticalLayout(fTabPage2));

    //fillings in tab0_page2
    TGVerticalFrame *fTabPage2Frame = new TGVerticalFrame(fTabPage2, 10, 10, kVerticalFrame);
    {
        fSettingText = new TGTextEdit(fTabPage2Frame, 800, 500);
        fSettingText->LoadBuffer("------> Ready <-------");
        fTabPage2Frame->AddFrame(fSettingText, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY, 2, 2, 2, 2));
    }
    fTabPage2->AddFrame(fTabPage2Frame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    // show mainwindow
    fTab0->SetTab(0);
    fTab0->Resize(fTab0->GetDefaultSize());
    this->AddFrame(fTab0, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    this->MapSubwindows();
    this->MapWindow();

    // set initial value
    if (_dataswitch)
    {
        fDataSwitchButton->SetOn();
        DataSwitch();
		//cout<<"there "<<endl;
    }
    fSettingText->LoadBuffer(gMyRootClass->GenerateSettingsText());
}

MyRootGui::~MyRootGui()
{
    // Destructor.
    CloseWindow();
}

void MyRootGui::CloseWindow()
{
    // Destructor.
    gMyRootClass->ReadSettings(fSettingText->GetText());
    gMyRootClass->SaveSettingsToEnv(env);
    delete gMyRootClass;
    DeleteWindow();
    gApplication->Terminate(0);
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
    if (fi.fFilename == NULL)
        return;
    fileDir = fi.fIniDir;
    filePath = fi.fFilename;
    fFileLabel1->SetText(filePath);
    env->SetValue("fileDir", fileDir);
    env->SetValue("filePath", filePath);
    env->SaveLevel(kEnvLocal);;
    cout << "----> Open " << fileDir << " " << filePath << endl;

    fAnaButton->SetEnabled(kTRUE);
    fHSlider1->SetEnabled(kFALSE);
    fPreButton->SetEnabled(kFALSE);
    fNextButton->SetEnabled(kFALSE);
    fDrawSFbutton->SetEnabled(kFALSE);

    fShowAllButton->SetEnabled(kFALSE);
    fShowIPButton->SetEnabled(kFALSE);
    fShowPolButton->SetEnabled(kFALSE);
}

void MyRootGui::OpenPedFile()
{
    static TString dir(".");
    TGFileInfo fi;
    fi.fFileTypes = filetypes;
    fi.fIniDir = StrDup(dir);
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (fi.fFilename == NULL)
        return;
    pedPath = fi.fFilename;
    fFileLabel2->SetText(pedPath);
    env->SetValue("pedPath", pedPath);
    env->SaveLevel(kEnvLocal);;

    //fNSigPedButton->SetEnabled(kTRUE);
    //fShowPedButton->SetEnabled(kFALSE);
}

//______________________________________________________________________________
//
void MyRootGui::Analysis()
{
    gMyRootClass->DataSwitch(fDataSwitchButton->IsOn());
    gMyRootClass->ReadSettings(fSettingText->GetText());
    gMyRootClass->SaveSettingsToEnv(env);
    fOutputText->LoadBuffer(gMyRootClass->GenerateSettingsOutput());
    
    env->SetValue("DataSwitch", fDataSwitchButton->IsOn());
    env->SaveLevel(kEnvLocal);;

    int retval;
    new TGMsgBox(gClient->GetRoot(), this,
                 "Select", "Analysis the whole directory (choose YES) or single file (choose NO):",
                 kMBIconQuestion, kMBYes | kMBNo, &retval);

    if (retval == 2)
        gMyRootClass->Analysis(filePath);
    else
        gMyRootClass->AnalysisDir(fileDir);


    fCanvas1->cd();
    fAnaButton->SetEnabled(kFALSE);
    fHSlider1->SetEnabled(kTRUE);
    fPreButton->SetEnabled(kTRUE);
    fNextButton->SetEnabled(kTRUE);
    fDrawSFbutton->SetEnabled(kTRUE);

    fShowAllButton->SetEnabled(kTRUE);
    fShowIPButton->SetEnabled(kTRUE);
    fShowPolButton->SetEnabled(kTRUE);

    fHSlider1->SetRange(0, gMyRootClass->GetNEvent());
    fHSlider1->SetPosition(gMyRootClass->GetIp());

    ShowAllButtonFunc();
    fTab0->Layout();
}

void MyRootGui::DrawPre()
{
    fCanvas1->cd();
    gMyRootClass->DrawPre("box");
    fCanvas1->Modified();
    fCanvas1->Update();
    fOutputText->LoadBuffer(gMyRootClass->GetInfo()->Data());
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
    fOutputText->LoadBuffer(gMyRootClass->GetInfo()->Data());
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
    fOutputText->LoadBuffer(gMyRootClass->GetInfo()->Data());

    fCanvas2->cd();
    gMyRootClass->DrawRaw();
    fCanvas2->Modified();
    fCanvas2->Update();
}

void MyRootGui::DrawSearchFrameFunc()
{
    fCanvas1->cd();
    gMyRootClass->DrawSearchFrameFunc();
    fCanvas1->Modified();
    fCanvas1->Update();
}

//______________________________________________________________________________
//
void MyRootGui::ShowAllButtonFunc()
{
    fCanvas1->cd();
    gMyRootClass->ShowAllButtonFunc(1);
    fCanvas1->Modified();
    fCanvas1->Update();

    fCanvas2->cd();
    gMyRootClass->ShowAllButtonFunc(2, "box");
    fCanvas2->Modified();
    fCanvas2->Update();
}

void MyRootGui::ShowIPButtonFunc()
{
    fCanvas1->cd();
    gMyRootClass->ShowIPButtonFunc(1);
    fCanvas1->Modified();
    fCanvas1->Update();

    fCanvas2->cd();
    gMyRootClass->ShowIPButtonFunc(2);
    fCanvas2->Modified();
    fCanvas2->Update();
}

void MyRootGui::ShowPolButtonFunc()
{
    fCanvas1->cd();
    gMyRootClass->ShowPolButtonFunc(1);
    fCanvas1->Modified();
    fCanvas1->Update();

    fCanvas2->cd();
    gMyRootClass->ShowPolButtonFunc(2);
    fCanvas2->Modified();
    fCanvas2->Update();
}

//______________________________________________________________________________

//______________________________________________________________________________
//
void MyRootGui::NSigPedButtonFunc()
{
    gMyRootClass->NSigPedButtonFunc(fNumberEntry31->GetNumber());

    fCanvas1->cd();
    gMyRootClass->AnalysisPed(pedPath);
    fCanvas1->Modified();
    fCanvas1->Update();

    fCanvas2->cd();
    gMyRootClass->ShowPedSigmaButtonFunc();
    fCanvas2->Modified();
    fCanvas2->Update();

    fNSigPedButton->SetEnabled(kFALSE);
    fShowPedButton->SetEnabled(kTRUE);
}

void MyRootGui::ShowPedButtonFunc()
{
    gMyRootClass->NSigPedButtonFunc(fNumberEntry31->GetNumber());

    fCanvas1->cd();
    gMyRootClass->ShowPedMeanButtonFunc();
    fCanvas1->Modified();
    fCanvas1->Update();

    fCanvas2->cd();
    gMyRootClass->ShowPedSigmaButtonFunc();
    fCanvas2->Modified();
    fCanvas2->Update();
}

//______________________________________________________________________________

//______________________________________________________________________________
//
void MyRootGui::DataSwitch()
{
    gMyRootClass->DataSwitch(fDataSwitchButton->IsOn());

    fNSigPedButton->SetEnabled(fDataSwitchButton->IsOn());
    fShowPedButton->SetEnabled(fDataSwitchButton->IsOn());
    fNumberEntry31->SetState(fDataSwitchButton->IsOn());
}

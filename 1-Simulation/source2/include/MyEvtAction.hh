//*********************************************
//  This is auto generated by G4gen 0.6
//                                  author:Qian

#ifndef MyEvtAction_h
#define MyEvtAction_h 1

#include "G4UserEventAction.hh"

class MyEvtAction : public G4UserEventAction
{
public:
  MyEvtAction();
  virtual ~MyEvtAction();

  virtual void BeginOfEventAction(const G4Event *);
  virtual void EndOfEventAction(const G4Event *);

private:
};

#endif
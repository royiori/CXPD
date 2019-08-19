#ifndef MyPhysicsList_h
#define MyPhysicsList_h 1

#include "G4VModularPhysicsList.hh"

class G4VPhysicsConstructor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MyPhysicsList : public G4VModularPhysicsList
{
public:
  MyPhysicsList();
  ~MyPhysicsList();

  virtual void ConstructParticle();
  virtual void ConstructProcess();
  void SetCuts();

protected:
  void AddParameterisation();

private:
  G4String fPsName;
  
  G4VPhysicsConstructor *fEMPhysicsList;
  G4VPhysicsConstructor *fEMOptPhysicsList;
  G4VPhysicsConstructor *fEMExtraPhysicsList;
  G4VPhysicsConstructor *fDecayPhysicsList;
  G4VPhysicsConstructor *fHadronElasticPhysicsList;
  G4VPhysicsConstructor *fHadronInelasticPhysicsList;
  G4VPhysicsConstructor *fStoppingPhysicsList;
  G4VPhysicsConstructor *fIonPhysicsList;
  G4VPhysicsConstructor *fNeutronPhysicsList;


};

#endif

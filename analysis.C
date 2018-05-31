#include "MyRootGui.h"

MyRootGui *gMyRootGui;

void analysis()
{
    gMyRootGui = new MyRootGui(gClient->GetRoot(), 1250, 580); 
}

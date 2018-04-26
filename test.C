#include "MyRootGui.h"

MyRootGui *gMyRootGui;

void test()
{
    gMyRootGui = new MyRootGui(gClient->GetRoot(), 850, 650); 
}

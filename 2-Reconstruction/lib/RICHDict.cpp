// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME libdIRICHDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "inc/MyRootGui.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *MyRootGui_Dictionary();
   static void MyRootGui_TClassManip(TClass*);
   static void delete_MyRootGui(void *p);
   static void deleteArray_MyRootGui(void *p);
   static void destruct_MyRootGui(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MyRootGui*)
   {
      ::MyRootGui *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::MyRootGui));
      static ::ROOT::TGenericClassInfo 
         instance("MyRootGui", "inc/MyRootGui.h", 27,
                  typeid(::MyRootGui), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &MyRootGui_Dictionary, isa_proxy, 0,
                  sizeof(::MyRootGui) );
      instance.SetDelete(&delete_MyRootGui);
      instance.SetDeleteArray(&deleteArray_MyRootGui);
      instance.SetDestructor(&destruct_MyRootGui);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MyRootGui*)
   {
      return GenerateInitInstanceLocal((::MyRootGui*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MyRootGui*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *MyRootGui_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::MyRootGui*)0x0)->GetClass();
      MyRootGui_TClassManip(theClass);
   return theClass;
   }

   static void MyRootGui_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_MyRootGui(void *p) {
      delete ((::MyRootGui*)p);
   }
   static void deleteArray_MyRootGui(void *p) {
      delete [] ((::MyRootGui*)p);
   }
   static void destruct_MyRootGui(void *p) {
      typedef ::MyRootGui current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MyRootGui

namespace {
  void TriggerDictionaryInitialization_RICHDict_Impl() {
    static const char* headers[] = {
"inc/MyRootGui.h",
0
    };
    static const char* includePaths[] = {
"/home/stuf/root6-build/include",
"/home/stuf/geant4.10.03.p03/examples/basic/geant4track/CXPD-master/2-Reconstruction/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RICHDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$inc/MyRootGui.h")))  MyRootGui;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RICHDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "inc/MyRootGui.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"MyRootGui", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RICHDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RICHDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RICHDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RICHDict() {
  TriggerDictionaryInitialization_RICHDict_Impl();
}

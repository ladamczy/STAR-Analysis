// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__starmIupc

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
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCRpsCluster.h"
#include "ReadFillPositionFile.h"
#include "MatchFillPosition.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_StUPCEvent(void *p = 0);
   static void *newArray_StUPCEvent(Long_t size, void *p);
   static void delete_StUPCEvent(void *p);
   static void deleteArray_StUPCEvent(void *p);
   static void destruct_StUPCEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StUPCEvent*)
   {
      ::StUPCEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StUPCEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StUPCEvent", ::StUPCEvent::Class_Version(), "StUPCEvent.h", 19,
                  typeid(::StUPCEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::StUPCEvent::Dictionary, isa_proxy, 4,
                  sizeof(::StUPCEvent) );
      instance.SetNew(&new_StUPCEvent);
      instance.SetNewArray(&newArray_StUPCEvent);
      instance.SetDelete(&delete_StUPCEvent);
      instance.SetDeleteArray(&deleteArray_StUPCEvent);
      instance.SetDestructor(&destruct_StUPCEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StUPCEvent*)
   {
      return GenerateInitInstanceLocal((::StUPCEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::StUPCEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_StUPCTrack(void *p = 0);
   static void *newArray_StUPCTrack(Long_t size, void *p);
   static void delete_StUPCTrack(void *p);
   static void deleteArray_StUPCTrack(void *p);
   static void destruct_StUPCTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StUPCTrack*)
   {
      ::StUPCTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StUPCTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StUPCTrack", ::StUPCTrack::Class_Version(), "StUPCTrack.h", 17,
                  typeid(::StUPCTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::StUPCTrack::Dictionary, isa_proxy, 4,
                  sizeof(::StUPCTrack) );
      instance.SetNew(&new_StUPCTrack);
      instance.SetNewArray(&newArray_StUPCTrack);
      instance.SetDelete(&delete_StUPCTrack);
      instance.SetDeleteArray(&deleteArray_StUPCTrack);
      instance.SetDestructor(&destruct_StUPCTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StUPCTrack*)
   {
      return GenerateInitInstanceLocal((::StUPCTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::StUPCTrack*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_StUPCBemcCluster(void *p = 0);
   static void *newArray_StUPCBemcCluster(Long_t size, void *p);
   static void delete_StUPCBemcCluster(void *p);
   static void deleteArray_StUPCBemcCluster(void *p);
   static void destruct_StUPCBemcCluster(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StUPCBemcCluster*)
   {
      ::StUPCBemcCluster *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StUPCBemcCluster >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StUPCBemcCluster", ::StUPCBemcCluster::Class_Version(), "StUPCBemcCluster.h", 11,
                  typeid(::StUPCBemcCluster), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::StUPCBemcCluster::Dictionary, isa_proxy, 4,
                  sizeof(::StUPCBemcCluster) );
      instance.SetNew(&new_StUPCBemcCluster);
      instance.SetNewArray(&newArray_StUPCBemcCluster);
      instance.SetDelete(&delete_StUPCBemcCluster);
      instance.SetDeleteArray(&deleteArray_StUPCBemcCluster);
      instance.SetDestructor(&destruct_StUPCBemcCluster);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StUPCBemcCluster*)
   {
      return GenerateInitInstanceLocal((::StUPCBemcCluster*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::StUPCBemcCluster*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_StUPCVertex(void *p = 0);
   static void *newArray_StUPCVertex(Long_t size, void *p);
   static void delete_StUPCVertex(void *p);
   static void deleteArray_StUPCVertex(void *p);
   static void destruct_StUPCVertex(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StUPCVertex*)
   {
      ::StUPCVertex *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StUPCVertex >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StUPCVertex", ::StUPCVertex::Class_Version(), "StUPCVertex.h", 11,
                  typeid(::StUPCVertex), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::StUPCVertex::Dictionary, isa_proxy, 4,
                  sizeof(::StUPCVertex) );
      instance.SetNew(&new_StUPCVertex);
      instance.SetNewArray(&newArray_StUPCVertex);
      instance.SetDelete(&delete_StUPCVertex);
      instance.SetDeleteArray(&deleteArray_StUPCVertex);
      instance.SetDestructor(&destruct_StUPCVertex);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StUPCVertex*)
   {
      return GenerateInitInstanceLocal((::StUPCVertex*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::StUPCVertex*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_StUPCTofHit(void *p = 0);
   static void *newArray_StUPCTofHit(Long_t size, void *p);
   static void delete_StUPCTofHit(void *p);
   static void deleteArray_StUPCTofHit(void *p);
   static void destruct_StUPCTofHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StUPCTofHit*)
   {
      ::StUPCTofHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StUPCTofHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StUPCTofHit", ::StUPCTofHit::Class_Version(), "StUPCTofHit.h", 12,
                  typeid(::StUPCTofHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::StUPCTofHit::Dictionary, isa_proxy, 4,
                  sizeof(::StUPCTofHit) );
      instance.SetNew(&new_StUPCTofHit);
      instance.SetNewArray(&newArray_StUPCTofHit);
      instance.SetDelete(&delete_StUPCTofHit);
      instance.SetDeleteArray(&deleteArray_StUPCTofHit);
      instance.SetDestructor(&destruct_StUPCTofHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StUPCTofHit*)
   {
      return GenerateInitInstanceLocal((::StUPCTofHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::StUPCTofHit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_StRPEvent(void *p = 0);
   static void *newArray_StRPEvent(Long_t size, void *p);
   static void delete_StRPEvent(void *p);
   static void deleteArray_StRPEvent(void *p);
   static void destruct_StRPEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StRPEvent*)
   {
      ::StRPEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StRPEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StRPEvent", ::StRPEvent::Class_Version(), "StRPEvent.h", 19,
                  typeid(::StRPEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::StRPEvent::Dictionary, isa_proxy, 4,
                  sizeof(::StRPEvent) );
      instance.SetNew(&new_StRPEvent);
      instance.SetNewArray(&newArray_StRPEvent);
      instance.SetDelete(&delete_StRPEvent);
      instance.SetDeleteArray(&deleteArray_StRPEvent);
      instance.SetDestructor(&destruct_StRPEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StRPEvent*)
   {
      return GenerateInitInstanceLocal((::StRPEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::StRPEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_StUPCRpsTrack(void *p = 0);
   static void *newArray_StUPCRpsTrack(Long_t size, void *p);
   static void delete_StUPCRpsTrack(void *p);
   static void deleteArray_StUPCRpsTrack(void *p);
   static void destruct_StUPCRpsTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StUPCRpsTrack*)
   {
      ::StUPCRpsTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StUPCRpsTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StUPCRpsTrack", ::StUPCRpsTrack::Class_Version(), "StUPCRpsTrack.h", 19,
                  typeid(::StUPCRpsTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::StUPCRpsTrack::Dictionary, isa_proxy, 4,
                  sizeof(::StUPCRpsTrack) );
      instance.SetNew(&new_StUPCRpsTrack);
      instance.SetNewArray(&newArray_StUPCRpsTrack);
      instance.SetDelete(&delete_StUPCRpsTrack);
      instance.SetDeleteArray(&deleteArray_StUPCRpsTrack);
      instance.SetDestructor(&destruct_StUPCRpsTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StUPCRpsTrack*)
   {
      return GenerateInitInstanceLocal((::StUPCRpsTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::StUPCRpsTrack*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_StUPCRpsTrackPoint(void *p = 0);
   static void *newArray_StUPCRpsTrackPoint(Long_t size, void *p);
   static void delete_StUPCRpsTrackPoint(void *p);
   static void deleteArray_StUPCRpsTrackPoint(void *p);
   static void destruct_StUPCRpsTrackPoint(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StUPCRpsTrackPoint*)
   {
      ::StUPCRpsTrackPoint *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StUPCRpsTrackPoint >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StUPCRpsTrackPoint", ::StUPCRpsTrackPoint::Class_Version(), "StUPCRpsTrackPoint.h", 14,
                  typeid(::StUPCRpsTrackPoint), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::StUPCRpsTrackPoint::Dictionary, isa_proxy, 4,
                  sizeof(::StUPCRpsTrackPoint) );
      instance.SetNew(&new_StUPCRpsTrackPoint);
      instance.SetNewArray(&newArray_StUPCRpsTrackPoint);
      instance.SetDelete(&delete_StUPCRpsTrackPoint);
      instance.SetDeleteArray(&deleteArray_StUPCRpsTrackPoint);
      instance.SetDestructor(&destruct_StUPCRpsTrackPoint);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StUPCRpsTrackPoint*)
   {
      return GenerateInitInstanceLocal((::StUPCRpsTrackPoint*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::StUPCRpsTrackPoint*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_StUPCRpsCluster(void *p = 0);
   static void *newArray_StUPCRpsCluster(Long_t size, void *p);
   static void delete_StUPCRpsCluster(void *p);
   static void deleteArray_StUPCRpsCluster(void *p);
   static void destruct_StUPCRpsCluster(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StUPCRpsCluster*)
   {
      ::StUPCRpsCluster *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StUPCRpsCluster >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StUPCRpsCluster", ::StUPCRpsCluster::Class_Version(), "StUPCRpsCluster.h", 14,
                  typeid(::StUPCRpsCluster), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::StUPCRpsCluster::Dictionary, isa_proxy, 4,
                  sizeof(::StUPCRpsCluster) );
      instance.SetNew(&new_StUPCRpsCluster);
      instance.SetNewArray(&newArray_StUPCRpsCluster);
      instance.SetDelete(&delete_StUPCRpsCluster);
      instance.SetDeleteArray(&deleteArray_StUPCRpsCluster);
      instance.SetDestructor(&destruct_StUPCRpsCluster);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StUPCRpsCluster*)
   {
      return GenerateInitInstanceLocal((::StUPCRpsCluster*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::StUPCRpsCluster*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr StUPCEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *StUPCEvent::Class_Name()
{
   return "StUPCEvent";
}

//______________________________________________________________________________
const char *StUPCEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StUPCEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *StUPCEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *StUPCEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr StUPCTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *StUPCTrack::Class_Name()
{
   return "StUPCTrack";
}

//______________________________________________________________________________
const char *StUPCTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StUPCTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *StUPCTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *StUPCTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr StUPCBemcCluster::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *StUPCBemcCluster::Class_Name()
{
   return "StUPCBemcCluster";
}

//______________________________________________________________________________
const char *StUPCBemcCluster::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCBemcCluster*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StUPCBemcCluster::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCBemcCluster*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *StUPCBemcCluster::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCBemcCluster*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *StUPCBemcCluster::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCBemcCluster*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr StUPCVertex::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *StUPCVertex::Class_Name()
{
   return "StUPCVertex";
}

//______________________________________________________________________________
const char *StUPCVertex::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCVertex*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StUPCVertex::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCVertex*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *StUPCVertex::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCVertex*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *StUPCVertex::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCVertex*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr StUPCTofHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *StUPCTofHit::Class_Name()
{
   return "StUPCTofHit";
}

//______________________________________________________________________________
const char *StUPCTofHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCTofHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StUPCTofHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCTofHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *StUPCTofHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCTofHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *StUPCTofHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCTofHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr StRPEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *StRPEvent::Class_Name()
{
   return "StRPEvent";
}

//______________________________________________________________________________
const char *StRPEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StRPEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StRPEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StRPEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *StRPEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StRPEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *StRPEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StRPEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr StUPCRpsTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *StUPCRpsTrack::Class_Name()
{
   return "StUPCRpsTrack";
}

//______________________________________________________________________________
const char *StUPCRpsTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StUPCRpsTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *StUPCRpsTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *StUPCRpsTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr StUPCRpsTrackPoint::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *StUPCRpsTrackPoint::Class_Name()
{
   return "StUPCRpsTrackPoint";
}

//______________________________________________________________________________
const char *StUPCRpsTrackPoint::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsTrackPoint*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StUPCRpsTrackPoint::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsTrackPoint*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *StUPCRpsTrackPoint::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsTrackPoint*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *StUPCRpsTrackPoint::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsTrackPoint*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr StUPCRpsCluster::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *StUPCRpsCluster::Class_Name()
{
   return "StUPCRpsCluster";
}

//______________________________________________________________________________
const char *StUPCRpsCluster::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsCluster*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StUPCRpsCluster::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsCluster*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *StUPCRpsCluster::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsCluster*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *StUPCRpsCluster::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StUPCRpsCluster*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void StUPCEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class StUPCEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StUPCEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(StUPCEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StUPCEvent(void *p) {
      return  p ? new(p) ::StUPCEvent : new ::StUPCEvent;
   }
   static void *newArray_StUPCEvent(Long_t nElements, void *p) {
      return p ? new(p) ::StUPCEvent[nElements] : new ::StUPCEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_StUPCEvent(void *p) {
      delete ((::StUPCEvent*)p);
   }
   static void deleteArray_StUPCEvent(void *p) {
      delete [] ((::StUPCEvent*)p);
   }
   static void destruct_StUPCEvent(void *p) {
      typedef ::StUPCEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::StUPCEvent

//______________________________________________________________________________
void StUPCTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class StUPCTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StUPCTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(StUPCTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StUPCTrack(void *p) {
      return  p ? new(p) ::StUPCTrack : new ::StUPCTrack;
   }
   static void *newArray_StUPCTrack(Long_t nElements, void *p) {
      return p ? new(p) ::StUPCTrack[nElements] : new ::StUPCTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_StUPCTrack(void *p) {
      delete ((::StUPCTrack*)p);
   }
   static void deleteArray_StUPCTrack(void *p) {
      delete [] ((::StUPCTrack*)p);
   }
   static void destruct_StUPCTrack(void *p) {
      typedef ::StUPCTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::StUPCTrack

//______________________________________________________________________________
void StUPCBemcCluster::Streamer(TBuffer &R__b)
{
   // Stream an object of class StUPCBemcCluster.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StUPCBemcCluster::Class(),this);
   } else {
      R__b.WriteClassBuffer(StUPCBemcCluster::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StUPCBemcCluster(void *p) {
      return  p ? new(p) ::StUPCBemcCluster : new ::StUPCBemcCluster;
   }
   static void *newArray_StUPCBemcCluster(Long_t nElements, void *p) {
      return p ? new(p) ::StUPCBemcCluster[nElements] : new ::StUPCBemcCluster[nElements];
   }
   // Wrapper around operator delete
   static void delete_StUPCBemcCluster(void *p) {
      delete ((::StUPCBemcCluster*)p);
   }
   static void deleteArray_StUPCBemcCluster(void *p) {
      delete [] ((::StUPCBemcCluster*)p);
   }
   static void destruct_StUPCBemcCluster(void *p) {
      typedef ::StUPCBemcCluster current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::StUPCBemcCluster

//______________________________________________________________________________
void StUPCVertex::Streamer(TBuffer &R__b)
{
   // Stream an object of class StUPCVertex.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StUPCVertex::Class(),this);
   } else {
      R__b.WriteClassBuffer(StUPCVertex::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StUPCVertex(void *p) {
      return  p ? new(p) ::StUPCVertex : new ::StUPCVertex;
   }
   static void *newArray_StUPCVertex(Long_t nElements, void *p) {
      return p ? new(p) ::StUPCVertex[nElements] : new ::StUPCVertex[nElements];
   }
   // Wrapper around operator delete
   static void delete_StUPCVertex(void *p) {
      delete ((::StUPCVertex*)p);
   }
   static void deleteArray_StUPCVertex(void *p) {
      delete [] ((::StUPCVertex*)p);
   }
   static void destruct_StUPCVertex(void *p) {
      typedef ::StUPCVertex current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::StUPCVertex

//______________________________________________________________________________
void StUPCTofHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class StUPCTofHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StUPCTofHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(StUPCTofHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StUPCTofHit(void *p) {
      return  p ? new(p) ::StUPCTofHit : new ::StUPCTofHit;
   }
   static void *newArray_StUPCTofHit(Long_t nElements, void *p) {
      return p ? new(p) ::StUPCTofHit[nElements] : new ::StUPCTofHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_StUPCTofHit(void *p) {
      delete ((::StUPCTofHit*)p);
   }
   static void deleteArray_StUPCTofHit(void *p) {
      delete [] ((::StUPCTofHit*)p);
   }
   static void destruct_StUPCTofHit(void *p) {
      typedef ::StUPCTofHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::StUPCTofHit

//______________________________________________________________________________
void StRPEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class StRPEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StRPEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(StRPEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StRPEvent(void *p) {
      return  p ? new(p) ::StRPEvent : new ::StRPEvent;
   }
   static void *newArray_StRPEvent(Long_t nElements, void *p) {
      return p ? new(p) ::StRPEvent[nElements] : new ::StRPEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_StRPEvent(void *p) {
      delete ((::StRPEvent*)p);
   }
   static void deleteArray_StRPEvent(void *p) {
      delete [] ((::StRPEvent*)p);
   }
   static void destruct_StRPEvent(void *p) {
      typedef ::StRPEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::StRPEvent

//______________________________________________________________________________
void StUPCRpsTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class StUPCRpsTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StUPCRpsTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(StUPCRpsTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StUPCRpsTrack(void *p) {
      return  p ? new(p) ::StUPCRpsTrack : new ::StUPCRpsTrack;
   }
   static void *newArray_StUPCRpsTrack(Long_t nElements, void *p) {
      return p ? new(p) ::StUPCRpsTrack[nElements] : new ::StUPCRpsTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_StUPCRpsTrack(void *p) {
      delete ((::StUPCRpsTrack*)p);
   }
   static void deleteArray_StUPCRpsTrack(void *p) {
      delete [] ((::StUPCRpsTrack*)p);
   }
   static void destruct_StUPCRpsTrack(void *p) {
      typedef ::StUPCRpsTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::StUPCRpsTrack

//______________________________________________________________________________
void StUPCRpsTrackPoint::Streamer(TBuffer &R__b)
{
   // Stream an object of class StUPCRpsTrackPoint.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StUPCRpsTrackPoint::Class(),this);
   } else {
      R__b.WriteClassBuffer(StUPCRpsTrackPoint::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StUPCRpsTrackPoint(void *p) {
      return  p ? new(p) ::StUPCRpsTrackPoint : new ::StUPCRpsTrackPoint;
   }
   static void *newArray_StUPCRpsTrackPoint(Long_t nElements, void *p) {
      return p ? new(p) ::StUPCRpsTrackPoint[nElements] : new ::StUPCRpsTrackPoint[nElements];
   }
   // Wrapper around operator delete
   static void delete_StUPCRpsTrackPoint(void *p) {
      delete ((::StUPCRpsTrackPoint*)p);
   }
   static void deleteArray_StUPCRpsTrackPoint(void *p) {
      delete [] ((::StUPCRpsTrackPoint*)p);
   }
   static void destruct_StUPCRpsTrackPoint(void *p) {
      typedef ::StUPCRpsTrackPoint current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::StUPCRpsTrackPoint

//______________________________________________________________________________
void StUPCRpsCluster::Streamer(TBuffer &R__b)
{
   // Stream an object of class StUPCRpsCluster.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StUPCRpsCluster::Class(),this);
   } else {
      R__b.WriteClassBuffer(StUPCRpsCluster::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StUPCRpsCluster(void *p) {
      return  p ? new(p) ::StUPCRpsCluster : new ::StUPCRpsCluster;
   }
   static void *newArray_StUPCRpsCluster(Long_t nElements, void *p) {
      return p ? new(p) ::StUPCRpsCluster[nElements] : new ::StUPCRpsCluster[nElements];
   }
   // Wrapper around operator delete
   static void delete_StUPCRpsCluster(void *p) {
      delete ((::StUPCRpsCluster*)p);
   }
   static void deleteArray_StUPCRpsCluster(void *p) {
      delete [] ((::StUPCRpsCluster*)p);
   }
   static void destruct_StUPCRpsCluster(void *p) {
      typedef ::StUPCRpsCluster current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::StUPCRpsCluster

namespace {
  void TriggerDictionaryInitialization_G__starmIupc_Impl() {
    static const char* headers[] = {
"StUPCEvent.h",
"StUPCTrack.h",
"StUPCBemcCluster.h",
"StUPCVertex.h",
"StUPCTofHit.h",
"StRPEvent.h",
"StUPCRpsTrack.h",
"StUPCRpsTrackPoint.h",
"StUPCRpsCluster.h",
"ReadFillPositionFile.h",
"MatchFillPosition.h",
0
    };
    static const char* includePaths[] = {
"/Users/leszekad/STAR-Analysis/test/include",
"/Applications/root_v6.16.00/include",
"/Users/leszekad/STAR-Analysis/test/build/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "G__starmIupc dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$StUPCEvent.h")))  StUPCEvent;
class __attribute__((annotate("$clingAutoload$StUPCTrack.h")))  StUPCTrack;
class __attribute__((annotate("$clingAutoload$StUPCBemcCluster.h")))  StUPCBemcCluster;
class __attribute__((annotate("$clingAutoload$StUPCVertex.h")))  StUPCVertex;
class __attribute__((annotate("$clingAutoload$StUPCTofHit.h")))  StUPCTofHit;
class __attribute__((annotate("$clingAutoload$StRPEvent.h")))  StRPEvent;
class __attribute__((annotate("$clingAutoload$StUPCRpsTrack.h")))  StUPCRpsTrack;
class __attribute__((annotate("$clingAutoload$StUPCRpsTrackPoint.h")))  StUPCRpsTrackPoint;
class __attribute__((annotate("$clingAutoload$StUPCRpsCluster.h")))  StUPCRpsCluster;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "G__starmIupc dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCRpsCluster.h"
#include "ReadFillPositionFile.h"
#include "MatchFillPosition.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"StRPEvent", payloadCode, "@",
"StUPCBemcCluster", payloadCode, "@",
"StUPCEvent", payloadCode, "@",
"StUPCRpsCluster", payloadCode, "@",
"StUPCRpsTrack", payloadCode, "@",
"StUPCRpsTrackPoint", payloadCode, "@",
"StUPCTofHit", payloadCode, "@",
"StUPCTrack", payloadCode, "@",
"StUPCVertex", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("G__star-upc",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_G__starmIupc_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_G__starmIupc_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_G__starmIupc() {
  TriggerDictionaryInitialization_G__starmIupc_Impl();
}

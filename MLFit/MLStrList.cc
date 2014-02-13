#include "MLFit/MLStrList.hh"
#include "TObjString.h"

ClassImp(MLStrList);

MLStrList::MLStrList(const char* s1, const char* s2, const char* s3, const char* s4, const char* s5, const char* s6, const char* s7, const char* s8) :
  TList()
{
  if (s1 != NULL) Add(new TObjString(s1));
  if (s2 != NULL) Add(new TObjString(s2));
  if (s3 != NULL) Add(new TObjString(s3));
  if (s4 != NULL) Add(new TObjString(s4));
  if (s5 != NULL) Add(new TObjString(s5));
  if (s6 != NULL) Add(new TObjString(s6));
  if (s7 != NULL) Add(new TObjString(s7));
  if (s8 != NULL) Add(new TObjString(s8));
}

MLStrList::MLStrList(TObject *obj1, TObject *obj2, TObject *obj3, TObject *obj4)
{
  if (obj1 != 0) Add(obj1);
  if (obj2 != 0) Add(obj2);
  if (obj3 != 0) Add(obj3);
  if (obj4 != 0) Add(obj4);
  
}

MLStrList::MLStrList(const MLStrList &strList) :
  TList()
{
  TIterator *iter = strList.MakeIterator();
  while (TObject *obj = iter->Next()) Add(obj);
  delete iter;
}

void MLStrList::Add(const char* s)
{
  TList::Add(new TObjString(s));
}

void MLStrList::Add(const MLStrList &strlist)
{
  TIterator *strIter = strlist.MakeIterator();
  while (TObjString *str = (TObjString*)strIter->Next()) {
    Add(str->GetName());
  }
}



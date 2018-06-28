#ifndef Auto_Rec
#define Auto_Rec
#pragma once
#include "TContour.h"
#include "TScene.h"
int AutoRecognition (TScene *SceneIn, char NameFileEtalon[],long *AllKmpObj, double stepxin, double stepyin,bool useKmpObj);
int Save_KML (char* sFileNameEtal, char* sFileNameKM, double stepxin, double stepyin);
#endif
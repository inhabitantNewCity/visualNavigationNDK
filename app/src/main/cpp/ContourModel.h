#ifndef ContourModel_HH
#define ContourModel_HH
#pragma once
#include "TScene.h"
// главная функция прослеживания контуров
void FindListContour (TShtrihPicture *shtrpic, TScene *SceneIn);
// прослеживание внешнего контура
void FindExternalContour (TShtrihPicture *shtrpic, TContour *extcontour, int NumStr, int NumShtr);
// прослеживание внутреннего контура
void FindInternalContour (TShtrihPicture *shtrpic, TContour *intercontour, int NumStr, int NumShtr);
// проверка на нахождение точки в списке контуров
bool CheckPointInListContour (TScene *wSceneIn, int x, int y);
#endif

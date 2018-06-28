#include <list>
#include "TShtrih.h"
//
// Created by Tmp on 28.06.2018.
//

#ifndef RECOGNATIONALGORITHM_RECOGNATIONFACED_H
#define RECOGNATIONALGORITHM_RECOGNATIONFACED_H
int BinarizationImage (unsigned char *pPixels, int width, int height, int threshold);
			 // перевод изображения в оттенки серого
void ImageToGrayscale (unsigned char *pPixels, size);
		 // определение порога методом Оцу
int OtsuThreshold(unsigned char *image, int size);
		 	// Функция создаёт список штрихов изображения
list<Shtrih> GoToLineFormat(unsigned char* mas, int width, int height);

int recognationTestTemplate(list<Shtrih> shtrihs, double stepxin, double stepyin,bool useKmpObj);
#endif //RECOGNATIONALGORITHM_RECOGNATIONFACED_H

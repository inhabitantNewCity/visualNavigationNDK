#include <jni.h>
#include <string>
#include "recognationFaced.h"

extern "C"
jstring
Java_navigation_pmk_com_recognationalgorithm_tetsActivity_stringFromJNI(
        JNIEnv* env,
        jobject /* this */) {
    std::string hello = "Hello from C++";
    return env->NewStringUTF(hello.c_str());
}

//please correct in parameters
jint Java_navigation_pmk_com_recognationalgorithm_tetsActivity_recognation(unsigned char *pPixels, int width, int height, int threshold){
    BinarizationImage (ppix, width, height, threshold);
    ImageToGrayscale(ppix, width);
    return recognationTestTemplate(GoToLineFormat(ppix, width, height),700, 700,false);
}
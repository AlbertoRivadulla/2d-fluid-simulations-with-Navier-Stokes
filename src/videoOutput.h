#ifndef VIDEO_OUTPUT_H
#define VIDEO_OUTPUT_H

#include <iostream>
#include <fstream>
#include <string>

#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>

// Create a video from the data in a file
void makeVideoFromDataFile( std::string inFileName, std::string outVideoName, int scale = 1 );

#endif

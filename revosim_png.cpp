#include <random>
#include <cmath>
#include<stdio.h>
#include <iostream>
#include "noise.h"
#include <typeinfo>
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>


using namespace std;

/* We want to be able to normalize noise from 0-255

   This function assumes a particular max (empiricaly obtained)
   and an assumed minimum of 0. If values exceed the new max,
   they are reset to 255.
   
   For pink noise the max is roughly 16
   
   MAX values are currently based on 10^6 values

*/

double maxmin_rgb(double m, double old_max){
		
	double new_max = 255; // sets the desired max
	double r;             // output
	
	r = (m/old_max)*new_max;  
	
	if(r > new_max){
		
		r = new_max;
		
	}
		
	return r;
}


/* BITMAP STRUCTS */

#define FHEADERSIZE 14
typedef struct{
	unsigned short type;
	unsigned int fSize;
	unsigned short reserved1;
	unsigned short reserved2;
	unsigned int offset;
}BMPFileHeader;

typedef struct{
	unsigned int hSize;
	int iWidth;
	int iHeight;
	unsigned short iPlanes;
	unsigned short iBPP; //Bits per pixel
	unsigned int iCompression;
	unsigned int iSize;
	int i_xPPM; //Pixel per Meter
	int i_yPPM;
	unsigned int iClrs;
	unsigned int iClrsImp;
}BMPInfoHeader;

typedef struct{
	unsigned char B;
	unsigned char G;
	unsigned char R;
}RGBColor; //24 Bit Color Struct



main(){
	
	/* 
	
	which_script == 0 calls for a text file, of n doubles, which is later passed to a python 
	script 'noisecheck.py'. That script fft's the noise and then plots the PSD. 
	
	which_script == 1, calculates the maximum amplitude of a signal of length n, of k generations
	of that signal. 
	
	which_script == 2, gets absolute of noise, then converts to 0-255 scale.
	
	*/
	
	int which_script = 4;                    //'bool' for functions
	signed long steps;
	std::cout << "Enter length of Signal: "; //pow(2,18);           //length of signal
	std::cin >> steps;
	double gam;
	std::cout << "Enter gamma: ";
	std::cin >> gam;						 //nature of signal
	std::cout << "*************************************" << "\n";
	std::vector<double> out;                 //allocating vector for signal
	
	if(which_script == 0){
	
		std::ofstream test;                
		
		test.open ("test.txt");

		out = noise(steps, gam, 1000, 0.01, 1, 3);

		for(signed long i = 0; i<steps; ++i){     //assigns each values from the vector to the txt file
			test << out[i] << " ";
		} 

		test.close();
	}
	
	if(which_script == 1){
		
		signed long num_signals = 10;
		
		double max = 0;
		
		//for(signed long i = 0; i<steps; i++){      // I used this to initially test absolute value was working.
			//std::cout << fabs(out[i]) << "\n";     // Seems to be working just fine. If not I think I can try,
		//}                                          // if val<0, -1*val, else val
		
		for(signed long i = 0; i<num_signals; i++){
			
			double temp = 0;
			out = noise(steps, gam, 1000, 0.01, 1, 3);
			
			for(signed long j = 0; j<steps; j++){

				if(fabs(out[j]) > temp){
					
					temp = fabs(out[j]);
					
				}		
			}
			
			std::cout << temp << "\n";
				
				
			if(temp > max){
				
				max = temp;
				
			}
		}
		
		std::cout << "Maximum Value: " << max << "\n";
		
		
	}
	
	if(which_script == 2){
		
		out = noise(32, gam, 1000, 0.01, 1, 3);
				
		for(signed long i = 0; i<32; i++){   //In short, this loop does abs(out[0:n-1])
			//std::cout << out[i] << "\n";
			out[i] = fabs(out[i]);
			std::cout << "Absolute Value: " << out[i] << "\n";
			std::cout << "0-255 Normalization: " << maxmin_rgb(out[i], gam) << "\n";
			//std::cout << out[i] << "\n";
		}
		
	}
		
		
		
		
		//std::cout << "TEST COMPLETE" << "\n";
		
	if(which_script == 3){
		
		signed long num_signals = 50;	//number of signals to be averaged
		
		std::vector<double> store_max; 
		store_max.reserve(num_signals);
		
		
		for(signed long k = 0; k<8; k++){ // (k+1)^10 == steps
			
			double avg = 0;
			
			double summation = 0;
			
			double sderr = 0;
			
			steps = steps*10;
			
			std::cout << steps << "\n";
		
			for(signed long i = 0; i<num_signals; i++){

				double temp = 0;
				out = noise(steps, gam, 1000, 0.01, 1, 3);

				for(signed long j = 0; j<steps; j++){

					if(fabs(out[j]) > temp){

						temp = fabs(out[j]);

					}		
				}

				avg += temp;
				store_max[i] = temp;

				//std::cout << temp << "\n";


				}
			
			avg = avg/num_signals;
			
			for(signed long i = 0; i<num_signals; i++){
				
				summation += pow(avg - store_max[i], 2);
				
			}
			
			sderr = (sqrt(summation/(num_signals - 1)))/sqrt(num_signals);

			std::cout << "Average value: " << avg << "\n";
			std::cout << "With Variance: " << sderr << "\n";
			
		}
		
	}
	
	if(which_script == 4){
		
		string pad;
		
		string filename;
		
		int r_val;
		int g_val;
		int b_val;
		
		double max;
		
		std::vector<double> out_r;
		std::vector<double> out_g;
		std::vector<double> out_b;
		
		out_r.reserve(steps);
		out_g.reserve(steps);
		out_b.reserve(steps);
		
		
		out_r = noise(steps, gam, 1000, 0.01, 1, 3);
		out_g = noise(steps, gam, 1000, 0.01, 1, 3);
		out_b = noise(steps, gam, 1000, 0.01, 1, 3);
		
		if(gam == 0){
			
			max = 1.96; //this is just, 1 sigma -- 97.5 percentile point
		}
		if(gam == 1){
		
			max = 1.6661*log10(steps) + 0.41152;
		}
		if(gam == 2){
		
			max = pow(10, 0.48078*log10(steps) + 0.14843);
		}
		
		for(signed long j = 0; j<steps; j++){
			
			if(j < 10){
				pad = "0000";
			} 
			if(j >= 10 && j < 100){
				pad = "000";
			}
			if(j >= 100 && j <1000){
				pad = "00";
			}
			if(j >= 1000 && j<10000){
				pad = "0";
			}
			if(j >= 10000){
				pad = "";
			}
			
			
			r_val = floor(maxmin_rgb(fabs(out_r[j]), max));
			g_val = floor(maxmin_rgb(fabs(out_g[j]), max));
			b_val = floor(maxmin_rgb(fabs(out_b[j]), max));
			
			
			string fext = ".bmp";

			ofstream bmpfile;
			filename = pad + to_string(j) + fext;
			bmpfile.open("/home/winston/noise_gen/Environment/" + filename, ios::out | ios::binary | ios::trunc);
			if(!bmpfile.is_open())
			{
				cout << "ERROR: FILE COULD NOT BE OPENED" << endl;
				return 1;
			}

			int width = 100;
			int height = 100;

			BMPFileHeader fHeader;   //File Header 
			fHeader.type = 19778;    //BM
			fHeader.fSize = FHEADERSIZE + sizeof(BMPInfoHeader) + sizeof(RGBColor)*width*height;   //File Size info
			fHeader.fSize = sizeof(fHeader) + sizeof(RGBColor);
			fHeader.reserved1 = 0;
			fHeader.reserved2 = 0;
			fHeader.offset = FHEADERSIZE + sizeof(BMPInfoHeader);    //Where we actually start writing the Bitmap data/image

			bmpfile.write((char*)(&fHeader.type), sizeof(fHeader.type));     //Let's start writing the file
			bmpfile.write((char*)(&fHeader.fSize), sizeof(fHeader.fSize));   //Header information
			bmpfile.write((char*)(&fHeader.reserved1), sizeof(fHeader.reserved1));
			bmpfile.write((char*)(&fHeader.reserved2), sizeof(fHeader.reserved2));
			bmpfile.write((char*)(&fHeader.offset), sizeof(fHeader.offset));

			BMPInfoHeader iHeader;
			iHeader.hSize = sizeof(BMPInfoHeader);
			iHeader.iWidth = width;
			iHeader.iHeight = height;
			iHeader.iPlanes = 1;
			iHeader.iBPP = 24;
			iHeader.iCompression = 0; //Uncompressed
			iHeader.iSize = width * height * 3;
			iHeader.i_xPPM = 0;
			iHeader.i_yPPM = 0;
			iHeader.iClrs = 0;
			iHeader.iClrsImp = 0;

			bmpfile.write((char*)&iHeader, sizeof(BMPInfoHeader));   //More Header information being placed in the file

			RGBColor* image;
			image = new RGBColor[width * height];
			
			for(int i = 0; i < width * height; i++) 
			{
				image[i].R = r_val;
				image[i].G = g_val;
				image[i].B = b_val;
			}
			//Write our arrays to the file
			for(int i = 0; i < width * height; i++)
			{
				bmpfile.write((char*)&image[i], sizeof(RGBColor));
			}
			// remove our temp image array

			delete [] image;
		   // close the file

			bmpfile.close();
			
			
		
		}
		
		return 0;
		
	}
	
}

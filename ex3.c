#include "Bitmap.h"
#define OUTPUT "sample_after.bmp"
#define ARYLEN(x) (sizeof(x)/sizeof(x[0]))

void printInfo(pBMPFILE p) {
  printBMPHeader(p);
  for (int i = 600; i < 600+10; i++) {
    printPixelAt(p, i);
  }
  return;
}

int main(int argc, char *argv[]) {
  pBMPFILE bmpfile;
  pBMPFILE output;
  pFILTER  filter;
  if(argc != 2) {
  perror("Usage: ex2 <target file>");
  return EXIT_FAILURE;
  }
  long gct[256];
  size_t pmax,pmin;
  bmpfile = prepareBMPFile(argv[1]);
  /*
  convertGrayScale(bmpfile);
  countPixel(gct,0,255,bmpfile);
  pmin = getPmin(gct, ARYLEN(gct), 50);
  pmax = getPmax(gct, ARYLEN(gct), 50);
  
  linearGreyLevelTransformation(bmpfile, (255.0/(pmax-pmin)), (-(255.0/(pmax-pmin)) * pmin));
  */
  output = copyBMPFile(bmpfile); 
  filter = newFilter();
  /*/
  setSlide(filter, LEFT_TO_RIGHT);
  for(int axj = 0; axj < 100; axj++) {
      if(axj%2) {
	  convolution(bmpfile,output,filter);
      }else{
	  convolution(output, bmpfile, filter);
      }
  }
  /*/
  setBokasi(filter);
  for(int axj = 0; axj < 10; axj++) {
      if(axj%2) {
	  convolution(bmpfile,output,filter);
      }else{
	  convolution(output, bmpfile, filter);
      }
  }// */
  writeBMPFile(output, OUTPUT);
  return 0;
}



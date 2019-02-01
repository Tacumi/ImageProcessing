#include "Bitmap.h"
#define OUTPUT "sample_after.bmp"

void printInfo(pBMPFILE p) {
  printBMPHeader(p);
  for (int i = 0; i < 10; i++) {
    printPixelAt(p, i);
  }
  return;
}

int main(int argc, char *argv[]) {
  pBMPFILE bmpfile;
  if(argc != 2) {
    perror("Usage: ex1 <target file>");
    return EXIT_FAILURE;
  }
  bmpfile = prepareBMPFile(argv[1]);
  // make yellow
  for(uint32_t i = 0; i < bmpfile->header.height; i++) {
    for(uint32_t j = 0; j < bmpfile->header.width; j++) {
      bmpfile->body[i*bmpfile->header.width + j].red = 0;
    }
  }
  
  printInfo(bmpfile);
  writeBMPFile(bmpfile, OUTPUT);
  return 0;
}



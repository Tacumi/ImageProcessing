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

uint8_t colorspace[256][256][256];

int main(int argc, char *argv[]) {
  pBMPFILE bmpfile;
  pBMPFILE tmp;
  pFILTER filter;
  pCLUSTER c;
  srand(time(NULL));
  if(argc != 2) {
    perror("Usage: ex4 <target file>");
    return EXIT_FAILURE;
  }
  c = malloc(CLUSTER_NUM*sizeof(CLUSTER));

  initCluster(c, CLUSTER_NUM);
 
  memset(colorspace, 0, sizeof(colorspace));
  bmpfile = prepareBMPFile(argv[1]);
  
 
  int f = 0;
  /*
  LOP(y, bmpfile->header.height) {
    LOP(x, bmpfile->header.width){
      PIXELMAP(bmpfile,x,y).red /= 128;
      PIXELMAP(bmpfile,x,y).red = PIXELMAP(bmpfile,x,y).red ? 191 : 63;
      PIXELMAP(bmpfile,x,y).blue /= 128;
      PIXELMAP(bmpfile,x,y).blue = PIXELMAP(bmpfile,x,y).blue ? 191 : 63;
      PIXELMAP(bmpfile,x,y).green /= 128;
      PIXELMAP(bmpfile,x,y).green = PIXELMAP(bmpfile,x,y).green ? 191 : 63;
    }
  }

  convertGrayScale(bmpfile);
  */
  for(int loop = 0; loop < 20; loop++){
    // assginment
    LOP(i,CLUSTER_NUM){
      c[i].pcount.red = 0;
      c[i].pcount.blue = 0;
      c[i].pcount.green = 0;
      c[i].count =0;
    }

    printf("\n");
    
    for(int y = 0; y < bmpfile->header.height; y++) {
      for(int x = 0; x < bmpfile->header.width; x++) {
	colorspace
	  [bmpfile->body->red]
	  [bmpfile->body->green]
	  [bmpfile->body->blue]
	  = getNearestClusterID(c,
				bmpfile->body[y*bmpfile->header.width + x],
				CLUSTER_NUM);
      }
    }
    LOP(i,CLUSTER_NUM){
      showCluster(c[i]);
    }
    // update
    LOP(idx, CLUSTER_NUM){
      if (c[idx].count != 0){
	c[idx].point.red = c[idx].pcount.red / c[idx].count;
	c[idx].point.blue = c[idx].pcount.blue / c[idx].count;
        c[idx].point.green = c[idx].pcount.green / c[idx].count;
      }
    }
  }
  mapColorSpace(c, bmpfile);
  tmp = copyBMPFile(bmpfile);
  /*
  filter = malloc(sizeof(FILTER));
  filter->value = malloc(sizeof(double*)*3);
  LOP(i,3){
    filter->value[i] = malloc(sizeof(double)*3);
  }

  setLaplacian4(filter);
  convolution(bmpfile,tmp,filter);
  
  LOP(y, bmpfile->header.height) {
    LOP(x, bmpfile->header.width){
      PIXELMAP(tmp,x,y).red = 255 - PIXELMAP(tmp,x,y).red;
      PIXELMAP(tmp,x,y).blue =  255 - PIXELMAP(tmp,x,y).blue;
      PIXELMAP(tmp,x,y).green =  255 - PIXELMAP(tmp,x,y).green;
    }
  }//*/
  writeBMPFile(tmp, OUTPUT);
  return 0;
}

#ifndef __BITMAP_H__
#define __BITMAP_H__

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define LOP(v,to) for(int v=0;v<(to);v++)
#define PIXELMAP(p,x,y) p->body[y*p->header.width + x]
#define CLUSTER_NUM 8

typedef struct _rgbspixel {
  uint8_t blue;
  uint8_t green;
  uint8_t red;
  uint8_t transparent;
} __attribute__((packed)) RGBSPIXEL, *pRGBSPIXEL;

typedef struct _rgbpixel {
  uint8_t blue;
  uint8_t green;
  uint8_t red;
} __attribute__((packed)) RGBPIXEL, *pRGBPIXEL;
typedef struct _filter{
  double **value;
  size_t width;
  size_t height;
} FILTER, *pFILTER;

typedef struct _graypixel {
  uint8_t gray;
} __attribute__((packed)) GRAYPIXEL, *pGRAYPIXEL;
struct rgbvalue{
  long long blue;
  long long green;
  long long red;
};
typedef struct _cluster {
  RGBPIXEL point;
  struct rgbvalue pcount;
  int32_t count;
  int32_t id;
} CLUSTER, *pCLUSTER;

typedef struct _bmpheader {
  char     identifier[2];
  uint32_t file_size;
  uint16_t reserve_0;
  uint16_t reserve_1;
  uint32_t header_size;
  uint32_t info_size;
  uint32_t width;
  uint32_t height;
  uint16_t num;
  uint16_t color_bit;
  uint32_t compress_algo;
  uint32_t compress_size;
  uint32_t hresolution;
  uint32_t vresolution;
  uint32_t color_num;
  uint32_t important_color_num;
} __attribute__((packed)) BMPHEADER, *pBMPHEADER;

typedef struct _bmpfile {
  BMPHEADER header;
  pRGBPIXEL   body;
} BMPFILE, *pBMPFILE;

enum DIRECTION{
	       LEFT_TO_RIGHT = 1,
	       RIGHT_TO_LEFT = 2,
	       UP_TO_DOWN    = 4,
	       DOWN_TO_UP    = 8,
};
// ex1
pBMPFILE newBMPBuffer();
pRGBPIXEL newBMPBody(uint32_t width, uint32_t height);
pBMPFILE readBMPHeader(pBMPFILE p, FILE *fp);
pBMPFILE readBMPBody(pBMPFILE p, FILE* fp);
int writeBMPFile(pBMPFILE p, char *write_to);
pBMPFILE prepareBMPFile(char *filename);
void printPixelAt(pBMPFILE p, int idx);
void printBMPHeader(pBMPFILE p);
// ex2
void countPixel(long *grayColorTable, int from, int to, pBMPFILE image);
void printGrayColorTable(long *grayColorTable, size_t size);
void convertGrayScale(pBMPFILE in);
size_t getPmin(long *table, size_t size, uint8_t threshold);
size_t getPmax(long *table, size_t size, uint8_t threshold);
pBMPFILE linearGreyLevelTransformation(pBMPFILE p, double a, double b);
uint8_t correctValue(double v);
// ex3
pBMPFILE convolution(pBMPFILE in, pBMPFILE out, pFILTER filter);
pBMPFILE copyBMPFile(pBMPFILE in);
pFILTER newFilter();
void setKrish(pFILTER filter, enum DIRECTION d);
void setLaplacian4(pFILTER p);
void setLaplacian8(pFILTER p);
void setInverseLaplacian4(pFILTER p);
void setInverseLaplacian8(pFILTER p);
void setthreshold(pBMPFILE p, uint8_t th);
void setSlide(pFILTER p, enum DIRECTION d);
void setBokasi(pFILTER p);
void setSobel(pFILTER p, enum DIRECTION d);
// ex4
double getDistanceRGB(RGBPIXEL a, RGBPIXEL b);
int32_t getNearestClusterID(pCLUSTER c, RGBPIXEL p, size_t siz);
void initCluster(pCLUSTER c, size_t siz);
void mapColorSpace(pCLUSTER c, pBMPFILE p);
void showCluster(CLUSTER c);
#endif /* __BITMAP_H__ */

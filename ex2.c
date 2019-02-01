#include "Bitmap.h"
#define OUTPUT "sample_after.bmp"
#define ARYLEN(x) (sizeof(x)/sizeof(x[0]))

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
    perror("Usage: ex2 <target file>");
    return EXIT_FAILURE;
  }
  // 画像を用意する．
  bmpfile = prepareBMPFile(argv[1]);
  long gct[256];
  size_t pmin, pmax;
  // ヒストグラムを用意する
  countPixel(gct,0,255, bmpfile);
  pmin = getPmin(gct, ARYLEN(gct), 0);
  pmax = getPmax(gct, ARYLEN(gct), 0);
  // 確認のために画面に表示
  printf("pmin : %ld\npmax : %ld\n",pmin,pmax);
  // ついでにヒストグラムを入れる
  printGrayColorTable(gct, ARYLEN(gct));
  //  convertGrayScale(bmpfile);
  double a = 255 / (pmax - (double)pmin);
  double b = -a * (double) pmin;
  // 線形変換を適用する
  linearGreyLevelTransformation(bmpfile, a, b);
  // もう一度カウントする．
  countPixel(gct,0,255, bmpfile);
  // ヒストグラムが広げられていたら0と255の位置になんらかの値が入っているはず．
  printGrayColorTable(gct, ARYLEN(gct));
  // ファイルに書き出す
  writeBMPFile(bmpfile, OUTPUT);
  return 0;
}



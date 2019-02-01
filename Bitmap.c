#include "Bitmap.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
uint32_t getImageRowSize(pBMPFILE p) {
  uint64_t ret;
  ret = p->header.color_bit;
  ret *= p->header.width;
  ret += 31;
  ret /= 32;
  ret *= 4;
  return ret;
}
pBMPFILE newBMPBuffer() {
  pBMPFILE ret = malloc(sizeof(BMPFILE) * 1);
  if (ret == NULL) {
    return NULL;
  }
  return ret;
}

pRGBPIXEL newBMPBody(uint32_t width, uint32_t height) {
  uint64_t siz = width * height;
  pRGBPIXEL ret = malloc(sizeof(RGBPIXEL) * siz);
  if (ret == NULL) {
    return NULL;
  }
  return ret;
}

pBMPFILE readBMPHeader(pBMPFILE p, FILE *fp) {
  fread(&(p->header), sizeof(BMPHEADER), 1, fp);
  return p;
}
/** readBMPBody(読み込ませるための構造体へのポインタ，ファイル名)
 *  画像本体を読み込む．パディングとか考えなきゃだめなのでfseekとか使ってずらす．
 */

pBMPFILE readBMPBody(pBMPFILE p, FILE* fp) {
  uint32_t imgx = p->header.width;
  uint32_t imgy = p->header.height;
  uint32_t padding = getImageRowSize(p) - sizeof(RGBPIXEL)*imgx;
  fseek(fp, 0, SEEK_SET);
  fseek(fp, p->header.header_size, SEEK_CUR);
  for(uint32_t y = 0; y < imgy; y++) {
    fread(p->body+(y*imgx), sizeof(RGBPIXEL), imgx, fp);
    fseek(fp, padding, SEEK_CUR);
  }
  return p;
}
/** writeBMPFile(読み込んだBMPファイルの構造体へのポインタ, 書き出す先のファイル名)
 *  ファイル名で指定したところに画像を書き出す
 */
int writeBMPFile(pBMPFILE p, char *write_to) {
  uint32_t imgx = p->header.width;
  uint32_t imgy = p->header.height;
  uint8_t EOL = 0;
  uint32_t padding = getImageRowSize(p) - sizeof(RGBPIXEL)*imgx;
  FILE *fp = fopen(write_to, "wb");
  fwrite(&(p->header), sizeof(BMPHEADER), 1, fp);
  for(uint32_t y = 0; y < imgy; y++) {
    fwrite(p->body+(y*imgx), sizeof(RGBPIXEL), imgx, fp);
    for(uint32_t i = 0; i < padding; i++) {
      fwrite(&EOL, sizeof(uint8_t), 1, fp);
    }
  }
  fclose(fp);
  return 0;
}
/** prepareBMPfile(ファイル名)
 *  ファイル名で指定したファイルから画像を読み込む
 */
pBMPFILE prepareBMPFile(char *filename){
  pBMPFILE bmpfile = newBMPBuffer();
  if (bmpfile == NULL) {
    exit(EXIT_FAILURE);
  }
  FILE* fp = fopen(filename, "rb");
  readBMPHeader(bmpfile,fp);
  bmpfile->body = newBMPBody(bmpfile->header.width, bmpfile->header.height);
  readBMPBody(bmpfile, fp);
  fclose(fp);
  return bmpfile;
}
/** printPixelAt(読み込んだ画像の構造体へのポインタ, 表示するインデックス(y*width+xで指定する))
 *  インデックスで指定したBMPファイルのピクセル情報を表示する
 */
void printPixelAt(pBMPFILE p, int idx) {
  printf("[Position %d]\n", idx);
  printf(" Blue:\t%2x", p->body[idx].blue);
  printf(" Green: %2x", p->body[idx].green);
  printf(" Red:\t%2x\n", p->body[idx].red);
  return;
}
/** countPixel(書き込むグレイスケール画像のヒストグラム，
 *             どの範囲からどの範囲までの色を指定するfrom,to，
 *             読み込んだ画像の構造体へのポインタ)
 *  グレイスケールにしてヒストグラムを取る
 *  RGBそれぞれに関してグレイスケールへの寄与度が違うのでそれに従って変換する
 */
void countPixel(long *grayColorTable, int from, int to, pBMPFILE image) {
  uint8_t gray;
  RGBPIXEL p;
  for (int i = 0; i < 255; i++) {
    grayColorTable[i] = 0;
  }
  for (uint32_t y = 0; y < image->header.height; y++) {
    for (uint32_t x = 0; x < image->header.width; x++) {
      p = image->body[y * image->header.width + x];
      gray = p.red * 0.3 + p.green * 0.59 + p.blue * 0.11;
      if(gray < from || gray > to){
        continue;
      }
      grayColorTable[gray] += 1;
	    
    }	    
  }
}
/** convertGrayScale(読み込んだBMPファイル構造体へのポインタ)
 *  BMPファイルをグレイスケールに変換する
 *  RGBそれぞれに関してグレイスケールへの寄与度が違うのでそれに従って変換する
 */
void convertGrayScale(pBMPFILE in) {
  uint8_t gray;
  pRGBPIXEL p;
  for (uint32_t y = 0; y < in->header.height; y++) {
    for (uint32_t x = 0; x < in->header.width; x++) {
      p = &(in->body[y * in->header.width + x]);
      gray = p->red * 0.3 + p->green * 0.59 + p->blue * 0.11;
      p->red = p->green = p->blue = gray;
    }	    
  }
}
/** printGrayColorTable(グレイスケールのヒストグラムを書き込むための配列，配列のサイズ)
 *  正しくヒストグラムが取れているかなーと確認するための関数
 */
void printGrayColorTable(long *grayColorTable, size_t size) {
  printf("GrayScale\n");
  for (size_t i = 0; i < size; i++) {
    printf("[%3ld]: %ld\n", i, grayColorTable[i]);
  }
}
/** printBMPHeader(読み込んだBMPファイル構造体へのポインタ)
 *  ヘッダ情報を表示するための関数．ファイルの読み込みとかは別の関数に任せる．
 *  TODO: NULLを渡されたときには落ちるのでなんとかする．(sitakunai)
 */
void printBMPHeader(pBMPFILE p) {
  pBMPHEADER ph = &(p->header);
  printf("Identifier: %c%c\n",
         ph->identifier[0],
         ph->identifier[1]);
  
  printf("File Size : %d\n",
         ph->file_size);

  printf("Reserved 0: %d\n",
         ph->reserve_0);

  printf("Reserved 1: %d\n",
         ph->reserve_1);

  printf("Header Size: %d\n",
         ph->header_size);

  printf("Info Size: %d\n",
         ph->info_size);

  printf("Width: %d\n",
         ph->width);

  printf("Height: %d\n",
         ph->height);

  printf("Num: %d\n",
         ph->num);

  printf("Colour Bit: %d\n",
         ph->color_bit);

  printf("Compress Algorithm: %d\n",
         ph->compress_algo);

  printf("Compressed Size: %d\n",
         ph->compress_size);

  printf("Horizontal Resolution: %d\n",
         ph->hresolution);

  printf("Vertical Resolution: %d\n",
         ph->vresolution);

  printf("Number of Colour: %d\n",
         ph->color_num);

  printf("Number of Important Colour: %d\n",
         ph->important_color_num);
  return;
}
/** getPmin(グレースケールに変換したときのヒストグラム, ヒストグラム配列の大きさ, 閾値)
    閾値以上のピクセル数を持つヒストグラムの最小値を求める
*/
size_t getPmin(long *table, size_t siz, uint8_t threshold) {
  for(int i = 0; i < siz; i++) {
    if(table[i] > threshold) {
      return i;
    }
  }
  return 255;
}
/** getPmax(グレースケールに変換したときのヒストグラム, ヒストグラム配列の大きさ, 閾値)
    閾値以上のピクセル数を持つヒストグラムの最大値を求める
*/
size_t getPmax(long *table, size_t siz, uint8_t threshold) {
  size_t ret = 0;
  for(int i = 0; i < siz; i++){
    if(table[i] > threshold) {
      ret = i;
    }
  }
  return ret;
}
/** linearGreyLevelTransformation(読み込んだBMP画像を利用する構造体へのポインタ, double a, double b)
 *  各レベルに対して線形変換を適用する関数．
 *  after = a(before) + bみたいな感じで変更してる．
 */
pBMPFILE linearGreyLevelTransformation(pBMPFILE p, double a, double b) {
  pRGBPIXEL pp = NULL;
  for(uint32_t y = 0; y < p->header.height; y++) {
    for(uint32_t x = 0; x < p->header.width; x++) {
      pp = (p->body + (y*p->header.width + x));
      pp->red = correctValue(pp->red * a + b);
      pp->green = correctValue(pp->green * a + b);
      pp->blue = correctValue(pp->blue * a + b);
    }
  }
  return p;
}
/** correctValue(計算後の画素値)
 *  入力された値を0から255の間にむりやり抑えたものを1バイトで返す
 *  そのまま計算したものを代入したらアンダーフローやオーバーフローが起きてしまったため実装．
 */
uint8_t correctValue(double v) {
  if(v < 0) { return -v>255?255:-v; }
  if(v > 255) { return 255; }
  return v;
}
/** newFilter()
 *  width,heightを受け取ってその大きさのフィルタの原型を返すつもりだったけど
 *  3x3しか利用しないので引数を省略する．
 */
pFILTER newFilter()
{
  pFILTER ret = (pFILTER)malloc(sizeof(FILTER));
  ret->width = 3;
  ret->height = 3;
  ret->value = (double**)malloc(sizeof(double*) * ret->height);
  for(uint32_t i = 0; i < ret->height; i++){
    ret->value[i] = (double*)malloc(sizeof(double) * ret->width);
    ret->value[i][0] = 0;
    ret->value[i][1] = 0;
    ret->value[i][2] = 0;
  }
  return ret;
}

pBMPFILE convolution(pBMPFILE in, pBMPFILE out, pFILTER filter) {
  int32_t x,y,i,j;
  for(y = 0; y < in->header.height; y++) {
    for(x = 0; x < in->header.width; x++) {
      int32_t red,blue,green;
      red = 0;
      blue = 0;
      green = 0;
      if(y>0 && y<in->header.height-1 && x>0 && x<in->header.width-1){
        for(i = -1; i <= 1; i++) {
          for(j = -1; j <= 1; j++) {
            red += (int32_t)in->body[((y+i)*in->header.width) + (x+j)].red
              * filter->value[i+1][j+1];
            blue += (int32_t)in->body[((y+i)*in->header.width) + (x+j)].blue
              * filter->value[i+1][j+1];
            green += (int32_t)in->body[((y+i)*in->header.width) + (x+j)].green
              * filter->value[i+1][j+1];
          }
        }
      }
      out->body[(y*out->header.width) + x].red = correctValue(red);
      out->body[(y*out->header.width) + x].blue = correctValue(blue);
      out->body[(y*out->header.width) + x].green = correctValue(green);
    }
  }
  return out;
}
pBMPFILE copyBMPFile(pBMPFILE in){
  pBMPFILE ret = newBMPBuffer();
  ret->header = in->header;
  ret->body = newBMPBody(in->header.width, in->header.height);
  for(size_t y = 0; y < in->header.height; y++) {
    for(size_t x = 0; x < in->header.width; x++) {
      ret->body[y*in->header.width + x] = in->body[y*in->header.width + x];
    }
  }
  return ret;
}
int deleteBMPBuffer(pBMPFILE in){
  if(in){
    if(in->body) {
      free(in->body);
    }
    free(in);
  }else{
    return -1;
  }
  return 0;
}
void setKrish(pFILTER filter, enum DIRECTION d) {
  for(int i=0; i < 3; i++) {
    for(int j=0; j < 3; j++) {
      if(d & LEFT_TO_RIGHT) {
        if(j == 0) {
          filter->value[i][j] -= 1;
        }else if(j == 2) {
          filter->value[i][j] += 1;
        }
      }
      if(d & UP_TO_DOWN) {
        if(i == 0) {
          filter->value[i][j] -= 1;
        }else if(i == 2) {
          filter->value[i][j] += 1;
        }
      }
      if(d & RIGHT_TO_LEFT) {
        if(j == 0) {
          filter->value[i][j] += 1;
        }else if(j == 2) {
          filter->value[i][j] -= 1;
        }

      }
      if(d & DOWN_TO_UP) {
        if(i == 0) {
          filter->value[i][j] += 1;
        }else if(i == 2) {
          filter->value[i][j] -= 1;
        }
      }
    }
  }
  for(int i=0; i < 3; i++) {
    for(int j=0; j < 3; j++) {
      if (filter->value[i][j] == 2 || filter->value[i][j] == -2){
        filter->value[i][j] /= 2;
      }
    }
  }
  return;
}
void setSobel(pFILTER p, enum DIRECTION d) {
  int32_t v[3][3]={{-1, 0, 1},
                   {-2, 0, 2},
                   {-1, 0, 1}};
  if(d == LEFT_TO_RIGHT){
    for(int y=0;y<3;y++){
      for(int x=0;x<3;x++){
        p->value[y][x]=v[y][x];
      }
    }
  }else{
    for(int y=0;y<3;y++){
      for(int x=0;x<3;x++){
        p->value[y][x]=v[x][y];
      }
    }
  }
}

void setLaplacian4(pFILTER p) {
  int32_t v[3][3]={{ 0,-1, 0},
                   {-1, 4,-1},
                   { 0,-1, 0}};
  for(int y=0;y<3;y++){
    for(int x=0;x<3;x++){
      p->value[y][x]=v[y][x];
    }
  }
}
void setInverseLaplacian4(pFILTER p) {
  int32_t v[3][3]={{ 0, 1, 0},
                   { 1,-4, 1},
                   { 0, 1, 0}};
  for(int y=0;y<3;y++){
    for(int x=0;x<3;x++){
      p->value[y][x]=v[y][x];
    }
  }
}
void setInverseLaplacian8(pFILTER p) {
  int32_t v[3][3]={{ 1, 1, 1},
                   { 1,-8, 1},
                   { 1, 1, 1}};
  for(int y=0;y<3;y++){
    for(int x=0;x<3;x++){
      p->value[y][x]=v[y][x];
    }
  }
}
void setLaplacian8(pFILTER p) {
  int32_t v[3][3]={{-1,-1,-1},
                   {-1, 8,-1},
                   {-1,-1,-1}};
  for(int y=0;y<3;y++){
    for(int x=0;x<3;x++){
      p->value[y][x]=v[y][x];
    }
  }
}
void setBokasi(pFILTER p){
  double v[3][3] = {{0.125, 0.125, 0.125},
                    {0.125,     0, 0.125},
                    {0.125, 0.125, 0.125}};
  for(int y=0;y<3;y++){
    for(int x=0;x<3;x++){
      p->value[y][x]=v[y][x];
    }
  }
  return;  
}
void setSlide(pFILTER p, enum DIRECTION d) {
  if(p->value[1][0] && p->value[0][1]) {
    p->value[0][0] = 1;
  }
  else if(d & UP_TO_DOWN & RIGHT_TO_LEFT) {
    p->value[0][2] = 1;
  }
  else if(d & DOWN_TO_UP & LEFT_TO_RIGHT) {
    p->value[2][0] = 1;
  }
  else if(d & DOWN_TO_UP & RIGHT_TO_LEFT) {
    p->value[2][2] = 1;
  }
  else if(d & LEFT_TO_RIGHT) {
    p->value[1][0] = 1;
  }
  else if(d & UP_TO_DOWN) {
    p->value[0][1] = 1;
  }
  else if(d & RIGHT_TO_LEFT) {
    p->value[1][2] = 1;
  }
  else if(d & DOWN_TO_UP) {
    p->value[2][1] = 1;
  }
  return;  
}
void setthreshold(pBMPFILE p, uint8_t th) {
  for(size_t y = 0; y < p->header.height; y++) {
    for(size_t x = 0; x < p->header.width; x++) {
      pRGBPIXEL pp = p->body;
      pp[y*p->header.width + x].red =
        pp[y*p->header.width + x].red > th ? 255 : 0;
      pp[y*p->header.width + x].blue =
        pp[y*p->header.width + x].blue > th ? 255 : 0;
      pp[y*p->header.width + x].green =
        pp[y*p->header.width + x].green > th ? 255 : 0;

    }
  }
}
double getDistanceRGB(RGBPIXEL a, RGBPIXEL b) {
  double distred = a.red - b.red;
  double distblue = a.blue - b.blue;
  double distgreen = a.green - b.green;
  distred = distred * distred;
  distblue = distblue * distblue;
  distgreen = distgreen * distgreen;
  return distred + distblue + distgreen;
}
int32_t getNearestClusterID(pCLUSTER c, RGBPIXEL p, size_t siz) {
  double min = 99999999999;
  int32_t mini = -1;
  double d;
  for(int32_t i = 0; i < (int32_t)siz; i++) {
    d = getDistanceRGB(c[i].point, p);
    if(d < min){
      min = d;
      mini = i;
    }
  }
  c[mini].pcount.red += p.red;
  c[mini].pcount.blue += p.blue;
  c[mini].pcount.green += p.green;
  c[mini].count++;
  return mini+1;
}

void initCluster(pCLUSTER c, size_t siz) {
  for(int i = 0; i < (int)siz; i++) {
    c[i].id = i+1;
    c[i].point.red   = rand() % 256;
    c[i].point.green = rand() % 256;
    c[i].point.blue  = rand() % 256;
    c[i].count = 0;
  }
  return;
}
void mapColorSpace(pCLUSTER c, pBMPFILE p) {
  LOP(x,p->header.width){
    LOP(y,p->header.height){
      uint8_t mapped = 
	(getNearestClusterID(c, PIXELMAP(p,x,y), CLUSTER_NUM)-1) * ((double)255/CLUSTER_NUM);
      PIXELMAP(p,x,y).red = PIXELMAP(p,x,y).blue = PIXELMAP(p,x,y).green = mapped;
    }
  }
}

void showCluster(CLUSTER c){
  printf("id: %d, count: %6d, R: %d, G: %d, B: %d\n", c.id, c.count, c.point.red, c.point.green, c.point.blue);
  return;
}

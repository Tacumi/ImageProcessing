#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>/* abs() */
#include <limits.h>/* INT_MAX */
#include <inttypes.h>/* uint8_t */
#define ATTR_MAX 32
struct ClusterNode
{
	float center;
	float old_center;
	int count;
	int sum;
};
struct Image_attribute
{
	int label;
	int clusterNumber;
	float area;
	float x_sum;
	float y_sum;
	float x_wt;
	float y_wt;
};
struct Image
{
	char magic[3];
	int width;
	int height;
	int max;
	int cmax;
	int cmin;
	uint8_t image[1000][1000];
	struct Image_attribute attr[ATTR_MAX];
};
int readImage(char *filename, struct Image *dst);
int writeImage(char *filename, struct Image *src);
int getLabelmaxarea(const struct Image *image, int rank);
int getDistance(struct Image *image, struct Image *template, int x, int y, int minD);
int getPatternCount(struct Image *img);
void searchTemplate(struct Image *image, struct Image *template);
void searchTemplate2(struct Image *image, struct Image *template, int *retx, int *rety);
void inverseImage(struct Image *target);
void contrastImage(struct Image *target);
void initImage(struct Image *target);
void initImageAttribute(struct Image_attribute *image_attribute);
void adjustImage(struct Image *, int alpha, int beta);
void setImageBlock(struct Image *imgout, int y, int x, uint8_t avg, int bsize);
void mozaicImage(struct Image *img, struct Image *imgout, int bsize);
void zoomImage(struct Image *img, struct Image *imgout, int bsize);
void zoomImageFlexible(struct Image *img, struct Image *imgout, double heightRate, double widthRate);
void smoothImage_avg(struct Image *img, struct Image *imgout);
void smoothImage_med(struct Image *img, struct Image *imgout);
void rotateImage(struct Image *img, struct Image *imgout, double theta, int cx, int cy);
void binarizationImage(struct Image *img, struct Image *imgout, uint8_t threshold);
void affinTransformImage(struct Image *img, struct Image *imgout, double a, double b, double c, double d, double e, double f);
void labelingImage(struct Image *img, struct Image *imgout);
void labelImageRecursive(struct Image *img, int i, int j, uint8_t label);
void analyzeImage(struct Image *img);
void printImageAttribute(const struct Image *img);
void printImageClusterNumber(const struct Image *img);
void applyMaskImage(struct Image *imgout, const struct Image *mask);
void repeatFunction(void(*func)(struct Image *img, struct Image *imgout), struct Image *img, struct Image *imgout, int times);
void inflateImage(struct Image *img, struct Image *imgout);
void deflateImage(struct Image *img, struct Image *imgout);
void drawRectangle(struct Image *imgout, int x, int y, int width, int height, uint8_t color);
void clusteringImage(struct Image *imgout, int clusterCount);
void recalculateColor(struct Image *image);
void *safeMalloc(size_t siz);
void showUsage(const char *argvzero);
uint8_t getavg(struct Image *img, int x, int y);
uint8_t getmed(struct Image *img, int x, int y);
uint8_t getavg_block(struct Image *img, int y, int x, int bsize);
uint8_t getLinearColor(struct Image *img, double y, double x);
uint8_t getThreshold(struct Image *img);
struct Image *getImageMaskByLabel(struct Image *imgout, int label);
struct Image *cropImage(struct Image *imgout, const struct Image *src, int x, int y, int width, int height);

int main(int argc, char* argv[])
{
	char *input = "default.pgm";
	char *inputtemplate = "carte.pgm";
	char *output = "output.pgm";
	char *outtemplate = "outputtmp.pgm";
	struct Image *image;
	struct Image *temporary;
	struct Image *templateimage;
	struct Image *imageout;
	struct Image *imageouttemplate;
	int x;
	int y;

	if (argc == 2)
	{
		input = argv[1];
	}
	else if (argc == 3)
	{
		input = argv[1];
		inputtemplate = argv[2];
	}
	else
	{
		showUsage(argv[0]);
		exit(EXIT_SUCCESS);
	}

	image = (struct Image*) malloc(sizeof(struct Image));
	temporary = (struct Image*) safeMalloc(sizeof(struct Image));
	templateimage = (struct Image*) safeMalloc(sizeof(struct Image));
	imageout = (struct Image*) safeMalloc(sizeof(struct Image));
	imageouttemplate = (struct Image*) safeMalloc(sizeof(struct Image));

	initImage(image);
	initImage(imageout);
	initImage(temporary);

	readImage(input,image);
	readImage(inputtemplate,templateimage);
	*temporary = *image;
	searchTemplate2(temporary,templateimage,&x,&y);
	writeImage(output,temporary);
	writeImage(outtemplate,cropImage(imageouttemplate,image,x,y,templateimage->width,templateimage->height));

	free(image);	
	free(temporary);
	free(templateimage);
	free(imageout);
	free(imageouttemplate);
	return 0;
}
void showUsage(const char *argvzero)
{
	fprintf(stderr,"Usage :\n");
	fprintf(stderr,"%s <input> [<template>]\n", argvzero);
	return;
}
struct Image *cropImage(struct Image *imgout, const struct Image *src, int x, int y, int width, int height)
{
	int i;
	int j;
	memcpy(imgout->magic,src->magic,sizeof(char)*strlen(src->magic)+1);
	imgout->max = src->max;
	imgout->width = width;
	imgout->height = height;
	for(i = 0; i+y < src->height; i++)
	{
		for(j = 0; j+x < src->width; j++)
		{
			imgout->image[i][j] = src->image[i+y][j+x];
		}
	}
	recalculateColor(imgout);
	return imgout;
}
void recalculateColor(struct Image *img)
{
	int i;
	int j;
	img->cmax = 0;
	img->cmin = 255;
	for (i = 0; i < img->height; i++)
	{
		for (j = 0; j < img->width; j++)
		{
			if (img->image[i][j] < img->cmin)
			{
				img->cmin = img->image[i][j];
			}
			if (img->image[i][j] > img->cmax)
			{
				img->cmax = img->image[i][j];
			}
		}
	}
	return;
}
void *safeMalloc(size_t siz)
{
	void *ret;
	ret = malloc(siz);
	if (ret == NULL)
	{
		fprintf(stderr, "can't allocate memory\n");
		exit(EXIT_FAILURE);
	}
	memset(ret, 0, siz);
	return ret;
}
int getPatternCount(struct Image *img)
{
	int i;
	for (i = 0; img->attr[i].area && i < ATTR_MAX; i++);
	return i;
}
void clusteringImage(struct Image *imgout, int clusterCount)
{
	int i;
	int p;
	int D;
	int sumdiff;
	int minD;
	const int patternCnt = getPatternCount(imgout);
	struct ClusterNode *clst = (struct ClusterNode*) safeMalloc(sizeof(struct ClusterNode) * clusterCount);
	for (i = 0; i < clusterCount; i++) clst[i].center = imgout->attr[i].area;
	do
	{
		sumdiff = 0;
		for (i = 0; i < clusterCount; i++)
		{
			clst[i].sum = 0;
			clst[i].count = 0;
		}
		for (p = 0; p < patternCnt; p++)
		{
			minD = INT_MAX;
			for (i = 0; i < clusterCount; i++)
			{
				D = abs(clst[i].center - imgout->attr[p].area);
				if (D < minD)
				{
					minD = D;
					imgout->attr[p].clusterNumber = i;
				}
			}
		}
		for (i = 0; i < patternCnt; i++)
		{
			clst[imgout->attr[i].clusterNumber].count += 1;
			clst[imgout->attr[i].clusterNumber].sum += imgout->attr[i].area;
		}
		for (i = 0; i < clusterCount; i++)
		{
			clst[i].old_center = clst[i].center;
			clst[i].center = clst[i].sum / clst[i].count;
		}
		for (i = 0; i < clusterCount; i++)
		{
			sumdiff += abs(clst[i].center - clst[i].old_center);
		}
	}
	while (sumdiff > 0);
	free(clst);
}
void searchTemplate2(struct Image *image, struct Image *template, int *retx, int *rety)
{
	double maxS = -1;
	double S;
	int matchx = 0;
	int matchy = 0;
	int i;
	int j;
	double IT = 0.0f;
	double II = 0.0f;
	double TT = 0.0f;
	int m;
	int n;
	for (m = 0; m < template->height; m++)
	{
		for (n = 0; n < template->width; n++)
		{
			TT += template->image[m][n] * (double) template->image[m][n];
		}
	}
	for (i = 0; i <= image->height - template->height; i++)
	{
		for (j = 0; j <= image->width - template->width; j++)
		{
			IT = 0.0f;
			II = 0.0f;
			for (m = 0; m < template->height; m++)
			{
				for (n = 0; n < template->width; n++)
				{
					IT += image->image[i + m][j + n] * (double) template->image[m][n];
					II += image->image[i + m][j + n] * (double) image->image[i + m][j + n];
				}
			}
			S = IT / (sqrt(II) * sqrt(TT));

			if (S > maxS)
			{
				maxS = S;
				matchx = j;
				matchy = i;
				printf("maxS = %f, matchx = %d, matchy = %d\n", maxS, matchx, matchy);
			}
		}
	}
	*retx = matchx;
	*rety = matchy;
	drawRectangle(image, matchx, matchy, template->width, template->height, template->max);
	return;
}
void searchTemplate(struct Image *image, struct Image *template)
{
	int minD = INT_MAX;
	int D;
	int matchx = 0;
	int matchy = 0;
	int i;
	int j;
	for (i = 0; i <= image->height - template->height; i++)
	{
		for (j = 0; j <= image->width - template->width; j++)
		{
			D = getDistance(image, template, j, i, minD);
			if (D < minD)
			{
				minD = D;
				matchx = j;
				matchy = i;
				printf("minD = %d, matchx = %d, matchy = %d\n", minD, matchx, matchy);
			}
		}
	}
	drawRectangle(image, matchx, matchy, template->width, template->height, template->max);
	return;
}
int getDistance(struct Image *image, struct Image *template, int x, int y, int minD)
{
	int ret = 0;
	int i;
	int j;

	for (i = 0; i < template->height; i++)
	{
		for (j = 0; j < template->width; j++)
		{
			ret += abs(image->image[y + i][x + j] - template->image[i][j]);
			if (minD < ret) return ret;
		}
	}
	return ret;
}
void drawRectangle(struct Image *imgout, int x, int y, int width, int height, uint8_t color)
{
	int i;
	for (i = y; i < y + height; i++)
	{
		imgout->image[i][x] = color;
		imgout->image[i][x + width] = color;

	}
	for (i = x; i < x + width; i++)
	{
		imgout->image[y][i] = color;
		imgout->image[y + height][i] = color;
	}
	return;
}

void repeatFunction(void(*func)(struct Image *input, struct Image *output), struct Image *img, struct Image *imgout, int times)
{
	int i;
	for (i = 0; i < times; i++)
	{
		(*func)(img, imgout);
		*img = *imgout;
	}
}
void applyMaskImage(struct Image *imgout, const struct Image *mask)
{
	int i;
	int j;
	for (i = 0; i < imgout->height; i++)
	{
		for (j = 0; j < imgout->width; j++)
		{
			imgout->image[i][j] *= mask->image[i][j];
		}
	}
	return;
}
struct Image *getImageMaskByLabel(struct Image *imgout, int label)
{
	int i;
	int j;
	struct Image *mask = malloc(sizeof(struct Image));
	*mask = *imgout;
	for (i = 0; i < imgout->height; i++)
	{
		for (j = 0; j < imgout->width; j++)
		{
			if (imgout->image[i][j] != label) mask->image[i][j] = 0;
			else mask->image[i][j] = 1;
		}
	}
	return mask;
}
int getLabelmaxarea(const struct Image *image, int rank)
{
	const struct Image_attribute *pattr[ATTR_MAX];
	const struct Image_attribute *tmp;
	int i;
	int j;
	int imax;
	int max;
	if (rank >= ATTR_MAX) return -1;
	for (i = 0; i < ATTR_MAX; i++)
	{
		pattr[i] = &(image->attr[i]);
	}
	for (i = 0; i < ATTR_MAX - 1; i++)
	{
		max = pattr[i]->area;
		imax = i;
		for (j = i + 1; j < ATTR_MAX; j++)
		{
			if (max < pattr[j]->area)
			{
				max = pattr[j]->area;
				imax = j;
			}
		}
		if (i != imax)
		{
			tmp = pattr[i];
			pattr[i] = pattr[imax];
			pattr[imax] = tmp;
		}
	}
	return pattr[rank]->label;
}
void deflateImage(struct Image *img, struct Image *imgout)
{
	int i, j;
	*imgout = *img;
	for (i = 1; i < imgout->height - 1; i++)
	{
		for (j = 1; j < imgout->width - 1; j++)
		{
			imgout->image[i][j] = img->image[i - 1][j - 1] & img->image[i - 1][j] & img->image[i - 1][j + 1] &
				img->image[i][j - 1]   &                      img->image[i][j + 1]   &
				img->image[i + 1][j - 1] & img->image[i + 1][j] & img->image[i + 1][j + 1] ;
		}
	}
}
void inflateImage(struct Image *img, struct Image *imgout)
{
	int i, j;
	*imgout = *img;
	for (i = 1; i < imgout->height - 1; i++)
	{
		for (j = 1; j < imgout->width - 1; j++)
		{
			imgout->image[i][j] = img->image[i - 1][j - 1] | img->image[i - 1][j] | img->image[i - 1][j + 1] |
				img->image[i][j - 1]   |                      img->image[i][j + 1]   |
				img->image[i + 1][j - 1] | img->image[i + 1][j] | img->image[i + 1][j + 1] ;
		}
	}
}
void printImageClusterNumber(const struct Image *img)
{
	int i;
	printf("label number\tarea\tcluster number\n");
	for (i = 0; i < ATTR_MAX; i++)
	{
		if (img->attr[i].area == 0) continue;
		printf("\t%d\t%.0f\t%d\n", i, img->attr[i].area, img->attr[i].clusterNumber);
	}
}
void printImageAttribute(const struct Image *img)
{
	int i;
	printf("label number\tarea\txwt\tywt\n");
	for (i = 0; i < ATTR_MAX; i++)
	{
		if (img->attr[i].area == 0) continue;
		printf("\t%d\t%.0f\t%.2f\t%.2f\n", i, img->attr[i].area, img->attr[i].x_wt, img->attr[i].y_wt);
	}
}
void analyzeImage(struct Image *img)
{
	int i;
	int j;
	for (i = 0; i < img->height; i++)
	{
		for (j = 0; j < img->width; j++)
		{
			img->attr[img->image[i][j]].area += 1;
			img->attr[img->image[i][j]].x_sum += j;
			img->attr[img->image[i][j]].y_sum += i;
		}
	}

	for (i = 0; i < ATTR_MAX; i++)
	{
		img->attr[i].x_wt = img->attr[i].x_sum / img->attr[i].area;
		img->attr[i].y_wt = img->attr[i].y_sum / img->attr[i].area;
	}
}
void initImageAttribute(struct Image_attribute *iat)
{
	iat->area = 0;
	iat->x_sum = 0;
	iat->y_sum = 0;
	iat->x_wt = 0;
	iat->y_wt = 0;
	iat->clusterNumber = 0;
}
void labelingImage(struct Image *img, struct Image *imgout)
{
	int label = 1;
	int i, j;
	*imgout = *img;
	for (i = 1; i < imgout->height - 1; i++)
	{
		for (j = 1; j < imgout->width - 1; j++)
		{
			if (imgout->image[i][j] == imgout->max)
			{
				labelImageRecursive(imgout, i, j, label);
				label++;
			}
		}
	}
	imgout->cmax = label - 1;
}
void labelImageRecursive(struct Image *img, int i, int j, uint8_t label)
{
	img->image[i][j] = label;
	if (i < 0 || i > img->height || j < 0 || j > img->width) return;
	if (img->image[i - 1][j - 1] == img->max)
		labelImageRecursive(img, i - 1, j - 1, label);
	if (img->image[i - 1][j] == img->max)
		labelImageRecursive(img, i - 1, j, label);
	if (img->image[i - 1][j + 1] == img->max)
		labelImageRecursive(img, i - 1, j + 1, label);
	if (img->image[i][j - 1] == img->max)
		labelImageRecursive(img, i, j - 1, label);
	if (img->image[i][j + 1] == img->max)
		labelImageRecursive(img, i, j + 1, label);
	if (img->image[i + 1][j - 1] == img->max)
		labelImageRecursive(img, i + 1, j - 1, label);
	if (img->image[i + 1][j] == img->max)
		labelImageRecursive(img, i + 1, j, label);
	if (img->image[i + 1][j + 1] == img->max)
		labelImageRecursive(img, i + 1, j + 1, label);
	return;
}
uint8_t getThreshold(struct Image *img)
{
	int *hist = (int*) malloc(sizeof(int) * img->max);
	float N = img->width * img->height;
	double omega = 0;
	double mew = 0;
	double mewt = 0;
	double max = 0.0f;
	uint8_t maxp = 0;
	double maximize;
	int i;
	int j;
	if (hist == NULL)
	{
		perror("malloc");
		exit(1);
	}
	memset(hist, 0, sizeof(int) *img->max);
	for (i = 0; i < img->height; i++)
	{
		for (j = 0; j < img->width; j++)
		{
			(*(hist + img->image[i][j]))++;
		}
	}
	for (i = 0; i < img->max; i++)
	{
		mewt += i * *(hist + i) / N;
	}
	for (i = 0; i < img->max; i++)
	{
		omega += *(hist + i) / N;
		mew += i * (*(hist + i)) / N;
		maximize = mewt * omega - mew;
		maximize *= maximize;
		maximize /= omega * (1 - omega);

		if (maximize > max)
		{
			max = maximize;
			maxp = i;
		}
	}
	free(hist);
	return maxp;
}
void binarizationImage(struct Image *img, struct Image *imgout, uint8_t threshold)
{
	int i;
	int j;

	*imgout = *img;
	initImage(imgout);
	for (i = 0; i < imgout->height; i++)
	{
		for (j = 0; j < imgout->width; j++)
		{
			imgout->image[i][j] = (img->image[i][j] > threshold) ? img->max : 0;
		}
	}
}
void affinTransformImage(struct Image *img, struct Image *imgout, double a, double b, double c, double d, double e, double f)
{
	int i;
	int j;
	uint8_t linerColor;
	double x0;
	double y0;
	const double multiplyer = 1 / (a * e - b * d);
	*imgout = *img;
	initImage(imgout);

	for (i = 0; i < imgout->height; i++)
	{
		for (j = 0; j < imgout->width; j++)
		{
			x0 = multiplyer * (j * e - b * i + b * f - c * e);
			y0 = multiplyer * (a * i - j * d + c * d - a * f);

			if (0 <= x0 && x0 < img->width - 1 &&
					0 <= y0 && y0 < img->height - 1)
			{
				linerColor = getLinearColor(img, y0, x0);
				imgout->image[i][j] = linerColor;
			}
		}
	}
}

uint8_t getLinearColor(struct Image *img, double y, double x)
{
	double ret;
	double alpha = x - (double)(int) x;
	double beta = y - (double)(int) y;
	int u = (int) x;
	int v = (int) y;

	ret = img->image[v][u] * (1 - alpha) * (1 - beta) +
		img->image[v][u + 1] * (alpha) * (1 - beta) +
		img->image[v + 1][u] * (1 - alpha) * (beta) +
		img->image[v + 1][u + 1] * (alpha) * (beta);

	return (uint8_t) ret;
}
void rotateImage(struct Image *img, struct Image *imgout, double rotate, int cx, int cy)
{
	int i;
	int j;
	double x0;
	double y0;
	uint8_t linerColor;
	const double c = cos(rotate * M_PI / 180.0);
	const double s = sin(rotate * M_PI / 180.0);
	*imgout = *img;
	initImage(imgout);

	for (i = 0; i < imgout->height; i++)
	{
		for (j = 0; j < imgout->width; j++)
		{
			x0 = (c * (j - cx)) + (s * (i - cy)) + cx;
			y0 = (( -s) * (j - cx)) + (c * (i - cy))  + cy;

			if (0 <= x0 && x0 < img->width - 1 &&
					0 <= y0 && y0 < img->height - 1)
			{
				linerColor = getLinearColor(img, y0, x0);
				imgout->image[i][j] = linerColor;
			}
		}
	}
}
void zoomImageFlexible(struct Image *img, struct Image *imgout, double heightRate, double widthRate)
{
	int i;
	int j;
	uint8_t linerColor;
	*imgout = *img;
	imgout->width = (imgout->width * widthRate);
	imgout->height = (imgout->height * heightRate);


	for (i = 0; i < imgout->height; i++)
	{
		for (j = 0; j < imgout->width; j++)
		{
			linerColor = getLinearColor(img, i / heightRate, j / widthRate);
			imgout->image[i][j] = linerColor;
		}
	}
}
void zoomImage(struct Image *img, struct Image *imgout, int bsize)
{
	int i;
	int j;
	uint8_t srcColor;
	*imgout = *img;
	imgout->width *= bsize;
	imgout->height *= bsize;

	for (i = 0; i < imgout->height; i += bsize)
	{
		for (j = 0; j < imgout->width; j += bsize)
		{
			srcColor = img->image[i / bsize][j / bsize];
			setImageBlock(imgout, i, j, srcColor, bsize);
		}
	}
}

void mozaicImage(struct Image *img, struct Image *imgout, int bsize)
{
	int i;
	int j;
	uint8_t avg;
	*imgout = *img;
	for (i = 0; i < img->height; i += bsize)
	{
		for (j = 0; j < img->width; j += bsize)
		{
			avg = getavg_block(img, i, j, bsize);
			setImageBlock(imgout, i, j, avg, bsize);
		}
	}
}
void setImageBlock(struct Image *imgout, int y, int x, uint8_t avg, int bsize)
{
	int i;
	int j;
	for (i = 0; i < bsize; i++)
	{
		for (j = 0; j < bsize; j++)
		{
			if (( y + i >= imgout->height) || (x + j >= imgout->width)) continue;
			imgout->image[y + i][x + j] = avg;
		}
	}
}

uint8_t getavg_block(struct Image *img, int y, int x, int bsize)
{
	int i;
	int j;
	int cnt = 0;
	int sum = 0;
	for (i = 0; i < bsize; i++)
	{
		for (j = 0; j < bsize; j++)
		{
			if (( y + i >= img->height) || (x + j >= img->width)) continue;
			sum += img->image[y + i][x + j];
			cnt++;
		}
	}
	return sum / cnt;
}
void smoothImage_med(struct Image *img, struct Image *imgout)
{
	int i;
	int j;
	*imgout = *img;
	for (i = 1; i < img->height - 1; i++)
	{
		for (j = 1; j < img->width - 1; j++)
		{
			imgout->image[i][j] = getmed(img, i, j);
		}
	}
	return;
}
void smoothImage_avg(struct Image *img, struct Image *imgout)
{
	int i;
	int j;
	*imgout = *img;
	for (i = 1; i < img->height - 1; i++)
	{
		for (j = 1; j < img->width - 1; j++)
		{
			imgout->image[i][j] = getavg(img, i, j);
		}
	}
	return;
}
uint8_t getmed(struct Image *img, int y, int x)
{
	uint8_t a[9];
	int idx = 0;
	int i;
	int j;
	uint8_t t;
	for (i = -1; i < 2; i++)
	{
		for (j = -1; j < 2; j++)
		{
			a[idx++] = img->image[y + i][x + j];
		}
	}
	for (i = 8; i > 0; i--)
	{
		for (j = 0; j < i; j++)
		{
			if (a[j] < a[j + 1])
			{
				t = a[j];
				a[j] = a[j + 1];
				a[j + 1] = t;
			}
		}
	}

	return a[4];
}
uint8_t getavg(struct Image *img, int y, int x)
{
	int sum = 0;
	int i;
	int j;
	for (i = -1; i < 2; i++)
	{
		for (j = -1; j < 2; j++)
		{
			sum += img->image[y + i][x + j];
		}
	}
	return sum / 9;
}

void contrastImage(struct Image *image)
{
	int i;
	int j;
	recalculateColor(image);
	for (i = 0; i < image->height; i++)
	{
		for (j = 0; j < image->width; j++)
		{
			image->image[i][j] = (( image->image[i][j] - image->cmin) * image->max) / (image->cmax - image->cmin);
		}
	}
}
void initImage(struct Image *image)
{
	int i;
	int j;
	for (i = 0; i < image->height; i++)
	{
		for (j = 0; j < image->width; j++)
		{
			image->image[i][j] = 0;
		}
	}
	for (i = 0; i < ATTR_MAX; i++)
	{
		image->attr[i].label = i;
		initImageAttribute(&(image->attr[i]));
	}
}

void inverseImage(struct Image *image)
{
	int i;
	int j;
	for (i = 0; i < image->height; i++)
	{
		for (j = 0; j < image->width; j++)
		{
			image->image[i][j] = image->max - image->image[i][j];
		}
	}
}
void adjustImage(struct Image *image, int alpha, int beta)
{
	int i;
	int j;
	for (i = 0; i < image->height; i++)
	{
		for (j = 0; j < image->width; j++)
		{
			image->image[i][j] = alpha * image->image[i][j] + beta;
			if (image->image[i][j] > 255) image->image[i][j] = 255;
		}
	}
}
int readImage(char *input, struct Image *image)
{
	FILE *infp;
	int i;
	int j;
	int32_t in;
	initImage(image);
	infp = fopen(input, "r");

	if (!infp)
	{
		fprintf(stderr, "can't open read file\n");
		return -1;
	}
	fscanf(infp, "%c%c", &image->magic[0], &image->magic[1]);
	image->magic[2] = 0;
	fscanf(infp, "%d", &image->width);
	fscanf(infp, "%d", &image->height);
	fscanf(infp, "%d", &image->max);
	image->cmax = 0;
	image->cmin = image->max;
	for (i = 0; i < image->height; i++)
	{
		for (j = 0; j < image->width; j++)
		{
			fscanf(infp, "%d", &in);
			image->image[i][j] = 0xff & in;
			if (image->image[i][j] > image->cmax) image->cmax = image->image[i][j];
			if (image->image[i][j] < image->cmin) image->cmin = image->image[i][j];
		}
	}
	fclose(infp);
	return 0;
}
int writeImage(char * output, struct Image *image)
{
	FILE *outfp;
	int i;
	int j;
	outfp = fopen(output, "w");
	if (!outfp)
	{
		fprintf(stderr, "can't open write file\n");
		return -1;
	}
	fprintf(outfp, "%s\n", image->magic);
	fprintf(outfp, "%d %d\n", image->width, image->height);
	fprintf(outfp, "%d\n", image->max);
	for (i = 0; i < image->height; i++)
	{
		for (j = 0; j < image->width; j++)
		{
			fprintf(outfp, "%d ", image->image[i][j]);
		}
	}
	fclose(outfp);
	return 0;
}

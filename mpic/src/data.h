#ifndef DATA_H
#define DATA_H
#define _DEFAULT_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <dirent.h> 
#include <ftw.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include <pcre.h>
#include <assert.h>
#include <time.h>
#include <linux/limits.h>		/* give MAX_PATH, which is unsafe !!!
					 * see -> http://insanecoding.blogspot.fr/2007/11/pathmax-simply-isnt.html
					 * I use it anyway : this code is not actually intended to be in production as it is
					 */
#include <malloc.h>

#include "tiffio.h"
#include "gdal.h"			/* used to copy projection system & geo coordinates from ref image to output image   					     */

#define MAX_FOLDER_DEPTH 1		/* limit the recursivity of folder exploration: ie: look for data in folders at MAX_FOLDER_DEPTH level under params.infolder */
#define OVECCOUNT 30 			/* should be a multiple of 3 -> allow to match OVECCOUNT/3 substrings with regexp (we need 2)                                */
#define NA_VAL 0			/* Not A Value val -> put 0 in it cause we do some adds                                                                      */
#define VERBOSE 3			/* Verbose level : control the number of messages printed by the programm, put 0 to turn off msgs                            */
// #define WINDOW_SIZE 11   		/* window size for the computational kernel:                                                                                 */
					/*          uncomment for faster computation (improve performance as it is known at compile time)			     */
					/*          comment to allow setting window size on the cmd line.                                                            */
					/* ! Don't forget to recompile after each code modification !                                                                */
#define BIN_TMP_FILE "outmat.dat"	/* temporary bin file where outmat chunks are written in.                                                                    */
#define MEAN_TMP_FILE "outmat_mean.dat"	/* temporary bin file where outmat_mean chunks are written in.                                                                    */
#define NACOUNT_TMP_FILE "outmat_nacount.dat"	/* temporary bin file where outmat_nacount chunks are written in.                                                                    */
#define VELOC_EW_TMP_FILE "outmat_veloc_ew.dat"	/* temporary bin file where outmat_nacount chunks are written in.                                                                    */
#define VELOC_NS_TMP_FILE "outmat_veloc_ns.dat"	/* temporary bin file where outmat_nacount chunks are written in.                                                                    */
// #define STATIC_DECOMP   		/* NBX & NBY are used to set the number of block along X & Y respectively. Uncomment to compute decomp at runtime	     */
#define SQUARE_DECOMP    		/* compute a domain decomposition with nb block along X = nb block along Y						     */
#define NBX 1
#define NBY 4
// #define CLEAR_IN_TASK_PRELOAD
#define STORE_KIJ // -> far less L1 cache misses
//#define STORE_IJK

#define FASTER_NACOUNT    		/* compute na_count_blck faster (should stay on, unless you want to modify this part of code ...)			     */
// #define FAST_HYPOT			/* use fast sqrt approximation inside hypot funct -> see in function code how to increase precision    			     */
// #define SECOND_NEWTON_PASS		/* do a second Newton's method pass to increase precision   								     */
// #define EVEN_FASTER			/* even faster sqrt computation : not very accurate, but speed-up by a factor 2. Usefull for quick testes or debug	     */

// value for cc=0 in micmac output
#define ZERO_CC_IN_TIF 128

#define NUM_FILES 4			/* number of files produced by micmac 4: Px1, Px2, Correl, Snr. Used for -l option					     */

/* SDATA:	macro helper to access the data in the image stack at col i, row j and layer k
 *		To browse the inner array (padding not included), loop on nrow, ncol, & nlayer
 *
 * SDATA_RAW:	macro helper to access the data in the stack at col i, row j and layer k, padding included (ie: (0,0,0) is in the padded zone)
 * 		To browse all the array (padding included), loop on height, width, & nlayer					 
 */
#ifdef STORE_IJK 
// data stored first along i, then along j and then along k
#define SDATA(s, i, j, k) s->data[(k) * s->stride_k + ((j) + s->row_padding) * s->width + ((i) + s->col_padding)]
#define SDATA_RAW(s, i, j, k) s->data[(k)*s->stride_k + (j)*s->width + (i)]
#endif
#ifdef STORE_KIJ				 
// data stored first along k, then along i and then along j
#define SDATA(s, i, j, k) s->data[((j) + s->row_padding) * s->stride_j + ((i) + s->col_padding) * s->nlayer + (k)]	
#define SDATA_RAW(s, i, j, k) s->data[(j)*s->stride_j + (i)*s->nlayer + (k)]
#endif

// macro helper to access the data in the image at col i and row j
 
// macro helper to access the data in the image at col i and row j, padding included (ie: (0,0,0) is in the padded zone)
//#define IDATA(s, i, j) s->data[(j)*s->width + (i)]
#define IDATA(s, i, j) s->data[((j) + s->row_padding) * s->width + ((i) + s->col_padding)]
#define IDATA_RAW(s, i, j) s->data[(j)*s->width + (i)]

// macro helper to test if (i, j) coordinates are in the image s bounds (image not padded)
#define INSIDE(s, i, j) (((i) >= 0) && ((i)  < (s->width)) && ((j) >= 0) && ((j)  < (s->height)))

// macro helpers for timings
#ifdef _OPENMP
#include "omp.h"
#define CLOCK_START(t) (t).start = omp_get_wtime();
#define CLOCK_END(t) (t).end = omp_get_wtime();
#define CLOCK_GET(t) ((float) ((t).end - (t).start))
#define CLOCK_SUM(t) (t).sum += ((float) ((t).end - (t).start));
#else
#define CLOCK_START(t) (t).start = clock()
#define CLOCK_END(t) (t).end = clock()
#define CLOCK_GET(t) (((float) ((t).end - (t).start)) / CLOCKS_PER_SEC)
#define CLOCK_SUM(t) (t).sum += ((float) ((t).end - (t).start)) / CLOCKS_PER_SEC
#endif
#define CLOCK_STOP(t) CLOCK_END(t); CLOCK_SUM(t)

struct time_mes {
#ifdef _OPENMP
	float start;
	float end;
#else
	clock_t start;
	clock_t end;
#endif
	float sum;
};

struct all_time {
	struct time_mes load;
	struct time_mes mean;
	struct time_mes process;
	struct time_mes write;
	struct time_mes clear;
	struct time_mes glob;
	struct time_mes loop;
	struct time_mes task_preload;
	struct time_mes task_process;
	struct time_mes preload_wait;
	struct time_mes process_wait;
};

struct stack_f {
	int nrow;
	int ncol;
	size_t nlayer;

	int col_padding;
	int row_padding;

	size_t width;
	size_t height;
	size_t stride_k;
	size_t stride_j;
	size_t stride_i;

	float *restrict data;
};

struct image_f {
	int nrow;
	int ncol;

	int col_padding;
	int row_padding;

	size_t width;
	size_t height;

	float *restrict data;
};

struct image_i {
	int nrow;
	int ncol;

	int col_padding;
	int row_padding;

	size_t width;
	size_t height;

	int *restrict data;
};

struct my_date {
	int day;
	int month;
	int year;
};

struct all_data {
	struct stack_f *restrict ew;
	struct stack_f *restrict ns;
	struct stack_f *restrict cc;

	struct image_i *restrict na_count;
	struct image_i *restrict na_count_blck;
	struct image_f *restrict mean_magn;
	struct image_f *restrict mean_veloc_ew;
	struct image_f *restrict mean_veloc_ns;
	struct image_f *restrict out_mat;

	struct my_date *dates_1;
	struct my_date *dates_2;

	float* duration;

	char** px1_paths;
	char** px2_paths;
	char** cc_paths;

	int* px1_is_tiled;
	int* px2_is_tiled;
	int* cc_is_tiled;
};

struct block_info {
	uint num_block_x;
	uint num_block_y;
	uint size_x;
	uint size_y;
	uint x_offset_start;
	uint x_offset_end;
	uint y_offset_start;
	uint y_offset_end;

	int first_x;
	int last_x;

	int first_y;
	int last_y;
};

struct domain_info {
	uint nb_block_x;
	uint nb_block_y;

	uint size_block_x;
	uint size_block_y;

	struct block_info *blocks;
};

struct parameters {
	int nlayer;
	int only_multitemp;
	int neighborhood;
	int output_mean;
	int output_nacount;
	int output_veloc;
	int max_num_layer;
	char *infolder;
	char *ref_image;
	char *output_dir;
	float available_GiB_mem;
	float xres;
	float yres;
	float px1_pas;
	float px2_pas;
	float min_displ;
	float min_cc;
	float min_matches;
};

int get_date(struct my_date* date, const char* date_str);
int same_date (const struct my_date* date1, const struct my_date* date2);

int allocate_domain(struct domain_info **domain, const int nb_block_x, const int nb_block_y, const int ncol, const int nrow, const int col_padding, const int row_padding);

struct all_data *allocate_all_data (const int nb_row, const int nb_col, const int nb_layer, const int column_padding, const int row_padding);
struct stack_f *allocate_stack_f (const int nb_row, const int nb_col, const int nb_layer, const int column_padding, const int row_padding);
struct image_f * allocate_image_f(const int nb_row, const int nb_col, const int column_padding, const int row_padding);
struct image_i * allocate_image_i(const int nb_row, const int nb_col, const int column_padding, const int row_padding);
void clear_image_i (struct image_i *image);
void clear_image_f (struct image_f *image);
void clear_stack_f (struct stack_f *stack);
void free_all_data (struct all_data *data);
void free_stack_f (struct stack_f *stack);
void free_image_f(struct image_f *image);
void free_image_i(struct image_i *image);
void free_domain(struct domain_info *domain);
#endif

#ifndef FUNC_H
#define FUNC_H

#include "data.h"

#define INT16_DT 1
#define FLOAT_DT 3
#define CC_DATA_TYPE uint8_t
#define OUT_DATA_TYPE uint8_t

int count_layers (const struct parameters *param);
int count_layers_impl (const int max_num_layer, const char *dirname, const int only_multitemp, const int level, const pcre *re);

int init_domain_decomposition(const struct parameters *param, struct domain_info **domain, const int avail_mem_gb, const int ncol, const int nrow, const int nlayer, const int col_padding, const int row_padding, int *output_in_file);
int compute_num_blocks(const struct parameters *param, const int nrow, const int ncol, const int nlayer,  const int col_padding, const int row_padding, const int avail_mem_gb, int *output_in_file);
unsigned long long int stack_size(const int nrow, const int ncol, const int nlayer, const int col_padding, const int row_padding);

int search_files_in_folder(struct all_data *data, const int nrow, const int ncol, const int padding, const struct parameters *param);
int search_files_impl(struct all_data *data, const int max_num_layer, const char *dirname, const char* ext, const int only_multitemp, const int level, int num_layer, const pcre *re);
int get_files_from_cmd_line(struct all_data *data, const int nlayer, const int argc, char *const *argv, const int opt_ind);

int load_stacks(struct all_data *restrict data, const struct domain_info *domain, const int ind_block_x, const int ind_block_y, const struct parameters *param, const int numthreads);
int load_px1(struct all_data *restrict data, const struct domain_info *domain, const int ind_block_x, const int ind_block_y, const float xres, const float px1_pas, const int numthreads);
int load_px2(struct all_data *restrict data, const struct domain_info *domain, const int ind_block_x, const int ind_block_y, const float yres, const float px2_pas, const int numthreads);
int load_cc(struct all_data *restrict data, const struct domain_info *domain, const int ind_block_x, const int ind_block_y, const int numthreads);
int count_na(struct all_data *restrict data, const float min_cc, const int numthreads);

int process_stacks(struct all_data *restrict data, const struct parameters *param, const int numthreads);
int compute_mean_magnitude(struct all_data *restrict data, const int numthreads);
int compute_mean_velocity(struct all_data *restrict data, const int numthreads);
void clear_stacks_first_touch(struct all_data *restrict data, const int numthreads);
void clear_stacks_padding_only(struct all_data *restrict data, const int numthreads);

int write_partial_results_in_bin_file(FILE *fileout, struct all_data *restrict data, struct domain_info *domain, const int i, const int j, const int nrow, const int ncol, const int mean_mat);
int write_partial_results(struct image_f *outmat, struct all_data *restrict data, struct domain_info *domain, const int i_blck, const int j_blck, const int numthreads, const int mean_mat);
int write_results_in_tiff_from_binfile(char* fileout, struct all_data *restrict data, struct domain_info *domain, const int ncol, const int nrow, const int mean_mat);
int write_results_in_tiff(char* fileout, struct image_f *outmat, struct all_data *restrict data, struct domain_info *domain);

int create_dir (const char *name);
int is_directory(const char* path);
int is_reg_file(const char* path);
int extract_dates(const char* dirname, struct my_date* date1, struct my_date* date2, const pcre *re);
const char *get_filename_ext(const char *filename);
char *get_full_name(char *fullname, const char *parent, const char *child);
void file_err_stop(const char* filename);
int get_image_info(const char* file, int* nrow, int* ncol, int* is_tiled);

float fast_hypotf(float x, float y);
#endif

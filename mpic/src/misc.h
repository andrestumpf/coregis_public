#ifndef MISC_H
#define MISC_H

#include "data.h"

int file_exists(char *file);
int mkdir_rec(char *path);
int get_cmd_line_options(struct parameters *params, const int argc, char *const *argv, int *path_opt_ind);
void print_usage(char *const *argv);
void print_help(char *const *argv);
void reset_timer(struct all_time *timings);
void print_time(struct all_time *timings);
void print_switches(void);
int copy_geotiff_metadata(const char* file_dest, const char* file_src);
#endif

#include "func.h"

unsigned long long int stack_size(const int nrow, const int ncol, const int nlayer, const int col_padding, const int row_padding)
{
	return sizeof(struct stack_f) + ((unsigned long long)nrow + 2 * row_padding) * ((unsigned long long)ncol + 2 * col_padding) * (unsigned long long)nlayer * sizeof(float);
}

int compute_num_blocks(const struct parameters *param, const int nrow, const int ncol, const int nlayer,  const int col_padding, const int row_padding, const int avail_mem_gb, int *output_in_file)
{

	unsigned long long int one_gb = 1024 * 1024 * 1024;
	unsigned long long int avail_mem = avail_mem_gb * one_gb;
	unsigned long long int mem_outmat = ncol * nrow * sizeof(float);
	if (param->output_mean) mem_outmat += mem_outmat;
	if (param->output_nacount) mem_outmat += mem_outmat;
	
	// compute memory needed to run in-core (all-in-memory, no tiling)
	unsigned long long int mem_needed_nocut = 0;
	unsigned long long int f_stack_fp = stack_size(nrow, ncol, nlayer, col_padding, row_padding);
	unsigned long long int f_image_fp = sizeof(struct image_f) + ((unsigned long long)nrow) * ((unsigned long long)ncol) * sizeof(float);
	unsigned long long int i_image_fp = sizeof(struct image_i) + ((unsigned long long)nrow) * ((unsigned long long)ncol) * sizeof(int);

	mem_needed_nocut += NUM_FILES * f_stack_fp; // EW NS CC
	mem_needed_nocut += 2 * i_image_fp; // NA.count + NA.count_blck
	mem_needed_nocut += 2 * f_image_fp; // mean.magn output.mat

	mem_needed_nocut += f_image_fp; // outmat_results
	if (param->output_mean) mem_needed_nocut += f_image_fp; // outmat_mean
	if (param->output_nacount) mem_needed_nocut += f_image_fp; // outmat_nacount

	mem_needed_nocut += NUM_FILES * nlayer * PATH_MAX * sizeof(char);	// paths
	mem_needed_nocut += 2 * nlayer * sizeof(struct my_date);	// dates

	if (VERBOSE) printf("Memory needed:\r\t\t\t%lli Go\n", mem_needed_nocut/(1024*1024*1024));
	
	// no need to cut -> 1 block
	if (avail_mem > mem_needed_nocut) {
		*output_in_file = 0;
		return 1;
	}

	// compute memory needed for in memory output matrix
	if (avail_mem < (mem_outmat + one_gb)) {
		*output_in_file = 1;
	} else {	// It remains at least 1GB memory for computing the blocks
		*output_in_file = 0;
		avail_mem -= mem_outmat;
	}

	// compute nb block for out-of-core processing
	unsigned long long int mem_used;
	int nb_block_along_dir = 0;
	do {
		nb_block_along_dir++;
#ifdef SQUARE_DECOMP
		int nrow_cut = (int) ceilf(nrow / (float)nb_block_along_dir);
		int ncol_cut = (int) ceilf(ncol / (float)nb_block_along_dir);
#else // cut along Y direction
		int nrow_cut = (int) ceilf(nrow / (float)nb_block_along_dir);
		int ncol_cut = ncol;
#endif
		mem_used = 	2 * NUM_FILES * stack_size(nrow_cut, ncol_cut, nlayer, col_padding, row_padding) // EW NS CC
				+ nrow_cut * ncol_cut * (2 * sizeof(float) + 2 * sizeof(int))	// mean, out, nacount, nacount_blck
				+ NUM_FILES * nlayer * PATH_MAX * sizeof(char)	// paths
				+ 2 * nlayer * sizeof(struct my_date);	// dates

	} while (mem_used > avail_mem);

	return nb_block_along_dir;
}

int init_domain_decomposition(const struct parameters *param, struct domain_info **domain, const int avail_mem_gb, const int ncol, const int nrow, const int nlayer, const int col_padding, const int row_padding, int *output_in_file)
{
	int nb_block_x, nb_block_y, min_num_blocks;

	min_num_blocks = compute_num_blocks(param, nrow, ncol, nlayer, col_padding, row_padding, avail_mem_gb, output_in_file);

#ifdef STATIC_DECOMP // use this at your own risk !
	nb_block_x = NBX;
	nb_block_y = NBY;
	*output_in_file = 0;
#else
	#ifdef SQUARE_DECOMP
		nb_block_x = min_num_blocks;
		nb_block_y = min_num_blocks;
	#else	// cut along Y direction
		nb_block_x = 1;
		nb_block_y = min_num_blocks;
	#endif
#endif
	return allocate_domain(domain, nb_block_x, nb_block_y, ncol, nrow, col_padding, row_padding);
}

void clear_stacks_first_touch(struct all_data *restrict data, const int numthreads)
{
#pragma omp parallel for schedule(runtime) shared(data) num_threads(numthreads)
	for (int j = 0; j< data->ew->height; j++) {
		for (int i = 0; i< data->ew->width; i++) {
			for (int k = 0; k< data->ew->nlayer; k++) {
				SDATA_RAW(data->ew, i, j, k) = 0.0;
				SDATA_RAW(data->ns, i, j, k) = 0.0;
				SDATA_RAW(data->cc, i, j, k) = 1.0;
				IDATA_RAW(data->na_count, i, j) = 0;
			}
			int i_img = i - data->ew->col_padding;
			int j_img = j - data->ew->row_padding;
			// if (INSIDE(data->na_count, i_img, j_img)) IDATA(data->na_count, i_img, j_img) = 0;
			
			if (INSIDE(data->mean_magn, i_img, j_img)) IDATA_RAW(data->mean_magn, i_img, j_img) = 0.0;
			if (INSIDE(data->out_mat, i_img, j_img)) IDATA_RAW(data->out_mat, i_img, j_img) = 0.0;
		}
	}
}

void clear_stacks_padding_only(struct all_data *restrict data, const int numthreads)
{
	int y_padding = data->ew->row_padding;
	int x_padding = data->ew->col_padding;
	#pragma omp parallel default(none) shared(data, x_padding, y_padding) num_threads(2)
	{
	#pragma omp single
	{
		#pragma omp task
		#pragma omp parallel for schedule(guided) default(none) shared(data, x_padding, y_padding) num_threads(numthreads - 1)
		for (int j = 0; j< data->ew->height; j++) {
			if (j < y_padding || j >= (data->ew->height - y_padding) ) { // up & down border
				for (int i = 0; i< data->ew->width; i++) {
					for (int k = 0; k< data->ew->nlayer; k++) {
						SDATA_RAW(data->ew, i, j, k) = 0.0;
						SDATA_RAW(data->ns, i, j, k) = 0.0;
						SDATA_RAW(data->cc, i, j, k) = 1.0; // this way computed na_count is 0 in paddings
					}
				}
			} else { // left & right border
				for (int i = 0; i < x_padding; i++) { // left border
					for (int k = 0; k< data->ew->nlayer; k++) {
						SDATA_RAW(data->ew, i, j, k) = 0.0;
						SDATA_RAW(data->ns, i, j, k) = 0.0;
						SDATA_RAW(data->cc, i, j, k) = 1.0;
					}
				}
				for (int i = (data->ew->width - x_padding); i < data->ew->width; i++) { // right border
					for (int k = 0; k< data->ew->nlayer; k++) {
						SDATA_RAW(data->ew, i, j, k) = 0.0;
						SDATA_RAW(data->ns, i, j, k) = 0.0;
						SDATA_RAW(data->cc, i, j, k) = 1.0;
					}
				}
			}
		}

		#pragma omp task
		clear_image_i(data->na_count);
		clear_image_i(data->na_count_blck);
	}
	}
}

int process_stacks(struct all_data *restrict data, const struct parameters *param, const int numthreads)
{
	int window_width;
	int window_height;
#ifndef WINDOW_SIZE
	window_width = data->ew->col_padding * 2 + 1;
	window_height = data->ew->row_padding * 2 + 1;
#else 	
	window_width = WINDOW_SIZE;
	window_height = WINDOW_SIZE;
#endif
	float max_na_blck = window_width * window_height * (float) data->ew->nlayer * param->min_matches;
	float max_na = ((float)data->ew->nlayer * (1.f - param->min_matches));


# pragma omp parallel for schedule(runtime) default(none) shared(data, param, window_width, window_height, max_na, max_na_blck) num_threads(numthreads)
	for (int j = 0; j< data->ew->nrow; j++) {
		for (int i = 0; i< data->ew->ncol; i++) {
			if ((IDATA(data->mean_magn, i, j) > param->min_displ) && (IDATA(data->na_count, i, j) <= max_na) && !(IDATA(data->na_count_blck, i, j) > max_na_blck)) {
				float sum_ew = 0.f;
				float sum_ns = 0.f;
				float sum_dist = 0.f;
				for (int k_blck = 0; k_blck < data->ew->nlayer; k_blck++) {
					for (int j_blck = 0; j_blck < window_height; j_blck++) {
						for (int i_blck = 0; i_blck < window_width; i_blck++) {
							sum_ew += SDATA_RAW(data->ew, i + i_blck, j + j_blck, k_blck);
							sum_ns += SDATA_RAW(data->ns, i + i_blck, j + j_blck, k_blck);
#ifdef FAST_HYPOT
							sum_dist += fast_hypotf(SDATA_RAW(data->ew, i + i_blck, j + j_blck, k_blck), SDATA_RAW(data->ns, i + i_blck, j + j_blck, k_blck));
#else
							sum_dist += hypotf(SDATA_RAW(data->ew, i + i_blck, j + j_blck, k_blck), SDATA_RAW(data->ns, i + i_blck, j + j_blck, k_blck));
#endif
						}
					}
				}
#ifdef FAST_HYPOT
				IDATA(data->out_mat, i, j) = fast_hypotf(sum_ew, sum_ns) / sum_dist;
#else
				IDATA(data->out_mat, i, j) = hypotf(sum_ew, sum_ns) / sum_dist;
#endif
			} else {
				IDATA(data->out_mat, i, j) = NA_VAL;
			}
		}
	}
	return 0;
}

int compute_mean_magnitude(struct all_data *restrict data, const int numthreads)
{
# pragma omp parallel for schedule(runtime) default(none) shared(data) num_threads(numthreads)
	for (int j = 0; j< data->ew->nrow; j++) {
		for (int i = 0; i< data->ew->ncol; i++) {
			float sum_ew = 0;
			float sum_ns = 0;
			for (int k = 0; k< data->ew->nlayer; k++) {
				sum_ew += SDATA(data->ew, i, j, k);
				sum_ns += SDATA(data->ns, i, j, k);
			}
			// here div by zero produces NaNs ...
			float mean_ew = sum_ew / (float) (data->ew->nlayer - IDATA(data->na_count, i, j));
			float mean_ns = sum_ns / (float) (data->ew->nlayer - IDATA(data->na_count, i, j));
#ifdef FAST_HYPOT
			IDATA(data->mean_magn, i, j) = fast_hypotf(mean_ew, mean_ns);
#else
			IDATA(data->mean_magn, i, j) = hypotf(mean_ew, mean_ns);
#endif
		}
	}
	return 0;
}

int compute_mean_velocity(struct all_data *restrict data, const int numthreads)
{
# pragma omp parallel for schedule(runtime) default(none) shared(data) num_threads(numthreads)
	for (int j = 0; j< data->ew->nrow; j++) {
		for (int i = 0; i< data->ew->ncol; i++) {
			float sum_ew = 0;
			float sum_ns = 0;
			for (int k = 0; k< data->ew->nlayer; k++) {
				sum_ew += SDATA(data->ew, i, j, k) / data->duration[k];
				sum_ns += SDATA(data->ns, i, j, k) / data->duration[k];
			}
			// here div by zero produces NaNs ...
			IDATA(data->mean_veloc_ew, i, j) = sum_ew / (float) (data->ew->nlayer - IDATA(data->na_count, i, j));
			IDATA(data->mean_veloc_ns, i, j) = sum_ns / (float) (data->ew->nlayer - IDATA(data->na_count, i, j));
		}
	}
	return 0;
}


int load_stacks(struct all_data *restrict data, const struct domain_info *domain, const int ind_block_x, const int ind_block_y, const struct parameters *param, const int numthreads)
{
	int error = 0;
	error += load_px1(data, domain, ind_block_x, ind_block_y, param->xres, param->px1_pas, numthreads);
	error += load_px2(data, domain, ind_block_x, ind_block_y, param->yres, param->px2_pas, numthreads);
	error += load_cc(data, domain, ind_block_x, ind_block_y, numthreads);
	error += count_na(data, param->min_cc, numthreads);
	return error;
}

int count_na(struct all_data *restrict data, const float min_cc, const int numthreads)
{
        for (int k = 0; k < data->ew->nlayer; k++) {
# pragma omp parallel for schedule(runtime) default(none) shared(k, data) num_threads(numthreads)
                for (int j = 0; j < data->ew->height; j++) {
                        for (int i = 0; i < data->ew->width; i++) {
                                if (SDATA_RAW(data->cc, i, j, k) < min_cc) {
	                        	IDATA_RAW(data->na_count, i, j)++;
	                                SDATA_RAW(data->ew, i, j, k) = (float) NA_VAL;
	                                SDATA_RAW(data->ns, i, j, k) = (float) NA_VAL;
                                }       
                        }       
                }       
        }
#ifndef FASTER_NACOUNT
# pragma omp parallel for schedule(runtime) default(none) shared(data) num_threads(numthreads)
        for (int j = 0; j < data->ew->nrow; j++) {
                for (int i = 0; i < data->ew->ncol; i++) {
        		int i_min = i - data->ew->col_padding;
			int i_max = i + data->ew->col_padding + 1;
			int j_min = j - data->ew->row_padding;
			int j_max = j + data->ew->row_padding + 1;
			int sum = 0;
			for (int j_blck = j_min; j_blck < j_max; j_blck++) {
				for (int i_blck = i_min; i_blck < i_max; i_blck++) {
					sum += IDATA(data->na_count, i_blck, j_blck);
				}
			}
			IDATA_RAW(data->na_count_blck, i, j) = sum;
                }       
        }
#else
        // do the same thing but faster: far less computation (number of operation approximately divided by a factor of window_size)
        // try to modify this code at your own risk ... it's quite ugly and tricky.
        // ...
        // if despite everything you really want to modify it:
        // the trick is here to compute the sum in the sliding window for the first j iteration of each thread
        // then as we move along i, just use the previous sum, remove the first column sum, and add the next column sum.
        // We do the same along j direction with rows.
        // We use col[] row[] arrays as rolling buffer for columns and rows sums.
        //
        // !!! Carefull with omp clause !!! if you change the scheduling (static, chunksize), it will not work anymore !
	int sum;
	int window_width = 2 * data->ew->col_padding + 1;
	int window_height = 2 * data->ew->row_padding + 1;
	int row[window_height];
	int col[window_width];
	int irow = 0;
	int icol = 0;
	int new_col, new_row, last_sum_i = 0, last_sum_j = 0;
	int sum_init = 0;
	int chunk_size = (int)ceilf((float)data->ew->nrow / (float)numthreads);

# pragma omp parallel for schedule(static,chunk_size) default(none) shared(chunk_size, window_height, window_width, data) firstprivate(sum_init, last_sum_i, last_sum_j) private(col, row, sum, new_col, new_row, irow, icol) num_threads(numthreads)
        for (int j = 0; j < data->ew->nrow; j++) {
        	if (j % chunk_size == 0) {
			int i_min = - data->ew->col_padding;
			int i_max = data->ew->col_padding + 1;
			int j_min = j - data->ew->row_padding;
			int j_max = j + data->ew->row_padding + 1;
			memset(col, 0, window_width * sizeof(int));
			memset(row, 0, window_height * sizeof(int));
			sum_init = 0;
			for (int j_blck = j_min; j_blck < j_max; j_blck++) {
				for (int i_blck = i_min; i_blck < i_max; i_blck++) {
					sum_init += IDATA(data->na_count, i_blck, j_blck);
					col[i_blck + data->ew->col_padding] += IDATA(data->na_count, i_blck, j_blck);
					row[j_blck - j + data->ew->row_padding] += IDATA(data->na_count, i_blck, j_blck);
				}
			}
		}
        	irow = (j + window_height - 1) % window_height;
                for (int i = 0; i < data->ew->ncol; i++) {
        		icol = (i + window_width - 1) % window_width;
	        	if (j % chunk_size == 0) {
	        		if (i == 0) {
	        			sum = sum_init;
	        			last_sum_i = sum;
	        			last_sum_j = sum;
	        		} else {
	        			sum = last_sum_i;
	        			sum -= col[icol]; // first col
	        			new_col = 0;
	        			for (int k = (j - data->ew->row_padding); k < (j + data->ew->row_padding + 1); k++)
	        				new_col += IDATA(data->na_count, i + data->ew->col_padding, k);
	        			col[icol] = new_col;
	        			sum += new_col;
	        			last_sum_i = sum;
	        		}
	        	} else {
	        		if (i == 0) {
	        			// recompute columns
	        			int i_min = i - data->ew->col_padding;
					int i_max = i + data->ew->col_padding + 1;
					int j_min = j - data->ew->row_padding;
					int j_max = j + data->ew->row_padding + 1;
					int icol_tmp = (icol + 1) % window_width;
					memset(col, 0, window_width * sizeof(int));
					for (int i_blck = i_min; i_blck < i_max; i_blck++) {
	        				for (int j_blck = j_min; j_blck < j_max; j_blck++) {
							col[icol_tmp] += IDATA(data->na_count, i_blck, j_blck);
						}
						icol_tmp = (icol_tmp + 1)%window_width;
					}
	        			sum = last_sum_j;
	        			sum -= row[irow]; // first row
	        			new_row = 0;
	        			for (int k = (i - data->ew->col_padding); k < (i + data->ew->col_padding + 1); k++)
	        				new_row += IDATA(data->na_count, k, j + data->ew->row_padding);
	        			row[irow] = new_row;
	        			sum += new_row;
	        			last_sum_i = sum;
	        			last_sum_j = sum;
	        		} else {
	        			sum = last_sum_i;
	        			sum -= col[icol]; // first col
	        			new_col = 0;
	        			for (int k = (j - data->ew->row_padding); k < (j + data->ew->row_padding + 1); k++)
	        				new_col += IDATA(data->na_count, i + data->ew->col_padding, k);
	        			col[icol] = new_col;
	        			sum += new_col;
	        			last_sum_i = sum;
	        		}
	        	}
			IDATA_RAW(data->na_count_blck, i, j) = sum;
                }
        }
#endif
        return 0;
}

int load_px1(struct all_data *restrict data, const struct domain_info *domain, const int ind_block_x, const int ind_block_y, const float xres, const float px1_pas, const int numthreads)
{
	int error = 0;
	int num_block = ind_block_y * domain->nb_block_x + ind_block_x;
#pragma omp parallel for schedule(runtime) default(none) shared(num_block, error, data, domain, stderr) num_threads(numthreads)
	for (int k = 0; k < data->ew->nlayer; k++) {
		if (data->px1_is_tiled[k]) {
			TIFF* tif_px1 = TIFFOpen(data->px1_paths[k], "r");
			if (tif_px1 != NULL) {
			        uint32 imagewidth, imagelength;
			        uint32 tilewidth, tilelength;
				uint16 config, datatype;
				uint32 i_tile, j_tile, i_buf, j_buf;
				tdata_t buf_px1;

				int i_init, j_init;				
				TIFFGetField(tif_px1, TIFFTAG_IMAGEWIDTH, &imagewidth);
				TIFFGetField(tif_px1, TIFFTAG_IMAGELENGTH, &imagelength);
				TIFFGetField(tif_px1, TIFFTAG_TILEWIDTH, &tilewidth);
				TIFFGetField(tif_px1, TIFFTAG_TILELENGTH, &tilelength);
				TIFFGetField(tif_px1, TIFFTAG_PLANARCONFIG, &config);
				TIFFGetField(tif_px1, TIFFTAG_DATATYPE, &datatype);
				// printf("PX1 type: %d -> %s\n",datatype, datatype==INT16_DT?"int 16":datatype==FLOAT_DT?"float 32":"unknown");
				buf_px1 = _TIFFmalloc(TIFFTileSize(tif_px1));
				// David : Quick hack to disable test on TIFFTAG_PLANARCONFIG tag. This is ok since we only use scalar data.
				if (config == PLANARCONFIG_CONTIG || 1) {
					int first_x = domain->blocks[num_block].first_x;
					int first_y = domain->blocks[num_block].first_y;
					int last_x = domain->blocks[num_block].last_x;
					int last_y = domain->blocks[num_block].last_y;
					long y_offset_start = domain->blocks[num_block].y_offset_start;
					long y_offset_end = domain->blocks[num_block].y_offset_end;
					j_init = 0;
					if (!first_y) {
						y_offset_start -= data->ew->row_padding;
						j_init -= data->ew->row_padding;
					}
					if (!last_y) {
						y_offset_end += data->ew->row_padding;
					}
					long tile_start_y = y_offset_start / tilelength;
					long tile_end_y = y_offset_end / tilelength;
					long offset_start_intile_y = y_offset_start - tile_start_y * tilelength;
					long offset_end_intile_y = y_offset_end - tile_end_y * tilelength;

					long x_offset_start = domain->blocks[num_block].x_offset_start;
					long x_offset_end = domain->blocks[num_block].x_offset_end;
					i_init = 0;
					if (!first_x) {
						x_offset_start -= data->ew->col_padding;
						i_init -= data->ew->col_padding;
					}
					if (!last_x) {
						x_offset_end += data->ew->col_padding;
					}
					long tile_start_x = x_offset_start / tilewidth;
					long tile_end_x = x_offset_end / tilewidth;
					long offset_start_intile_x = x_offset_start - tile_start_x * tilewidth;
					long offset_end_intile_x = x_offset_end - tile_end_x * tilewidth;
					// assumption : image is not compressed : random access is not possible for compressed tiff
					// todo : test compression
					long shift_tile_y = 0;
					for (j_tile = tile_start_y; j_tile <= tile_end_y; j_tile++) {
						long offset_start_buf_y = (j_tile == tile_start_y) ? offset_start_intile_y : 0;
						long offset_end_buf_y = (j_tile == tile_end_y) ? offset_end_intile_y : tilelength;

						long shift_tile_x = 0;
						for (i_tile = tile_start_x; i_tile <= tile_end_x; i_tile++) {
							long offset_start_buf_x = (i_tile == tile_start_x) ? offset_start_intile_x : 0;
							long offset_end_buf_x = (i_tile == tile_end_x) ? offset_end_intile_x : tilewidth;
							TIFFReadTile(tif_px1, buf_px1, i_tile * tilewidth, j_tile * tilelength, 0, 0);
							int inbuf_j = 0;
							switch(datatype) {
							case INT16_DT:
							{
								int16_t *px1_p = buf_px1;
								for (j_buf = offset_start_buf_y; j_buf < offset_end_buf_y; j_buf++) {
									int inbuf_i = 0;
									for (i_buf = offset_start_buf_x; i_buf < offset_end_buf_x; i_buf++) {
										int i = i_init + shift_tile_x + inbuf_i;
										int j = j_init + shift_tile_y + inbuf_j;
										SDATA(data->ew, i, j, k) = (float) (px1_p[j_buf * tilewidth + i_buf] * xres * px1_pas);
										inbuf_i++;
									}
									inbuf_j++;
								}
								break;
							}
							case FLOAT_DT:
							{
								float *px1_p = buf_px1;
								for (j_buf = offset_start_buf_y; j_buf < offset_end_buf_y; j_buf++) {
									int inbuf_i = 0;
									for (i_buf = offset_start_buf_x; i_buf < offset_end_buf_x; i_buf++) {
										int i = i_init + shift_tile_x + inbuf_i;
										int j = j_init + shift_tile_y + inbuf_j;
										SDATA(data->ew, i, j, k) = (float) (px1_p[j_buf * tilewidth + i_buf] * xres * px1_pas);
										inbuf_i++;
									}
									inbuf_j++;
								}
								break;
							}
							default:
								error++;
								break;
							}
							shift_tile_x += (tilewidth - offset_start_buf_x);
						}
						shift_tile_y += (tilelength - offset_start_buf_y);
					}
				} else if (config == PLANARCONFIG_SEPARATE) {
					fprintf(stderr, "px1 file is layered: case not coded yet\n");
					error++;
				}
				_TIFFfree(buf_px1);
				TIFFClose(tif_px1);
			} else {
				fprintf(stderr, "error opening px1 file\n");
				error++;
			}
		} else {
			TIFF* tif_px1 = TIFFOpen(data->px1_paths[k], "r");
			if (tif_px1 != NULL) {
				uint32 imagelength;
				uint16 config, datatype;
				uint32 ii, jj;
				tdata_t buf_px1;

				int i, j;				
				TIFFGetField(tif_px1, TIFFTAG_IMAGELENGTH, &imagelength);
				TIFFGetField(tif_px1, TIFFTAG_PLANARCONFIG, &config);
				TIFFGetField(tif_px1, TIFFTAG_DATATYPE, &datatype);
				// printf("PX1 type: %d -> %s\n",datatype, datatype==INT16_DT?"int 16":datatype==FLOAT_DT?"float 32":"unknown");

				buf_px1 = _TIFFmalloc(TIFFScanlineSize(tif_px1));
				if (config == PLANARCONFIG_CONTIG || 1) {
					int first_x = domain->blocks[num_block].first_x;
					int first_y = domain->blocks[num_block].first_y;
					int last_x = domain->blocks[num_block].last_x;
					int last_y = domain->blocks[num_block].last_y;
					int y_offset_start = domain->blocks[num_block].y_offset_start;
					int y_offset_end = domain->blocks[num_block].y_offset_end;
					j = 0;
					if (!first_y) {
						y_offset_start -= data->ew->row_padding;
						j -= data->ew->row_padding;
					}
					if (!last_y) {
						y_offset_end += data->ew->row_padding;
					}
					assert(y_offset_end <= imagelength);
					// assumption : image is not compressed : random access is not possible for compressed tiff
					// todo : test compression
					for (jj = y_offset_start; jj < y_offset_end; jj++) {
						TIFFReadScanline(tif_px1, buf_px1, jj, 0);
						int x_offset_start = domain->blocks[num_block].x_offset_start;
						int x_offset_end = domain->blocks[num_block].x_offset_end;
						i = 0;
						if (!first_x) {
							x_offset_start -= data->ew->col_padding;
							i -= data->ew->col_padding;
						}
						if (!last_x) {
							x_offset_end += data->ew->col_padding;
						}
						switch(datatype) {
						case INT16_DT:
						{
							int16_t *px1_p = buf_px1;
							for (ii = x_offset_start; ii < x_offset_end; ii++) {
								SDATA(data->ew, i, j, k) = (float) (px1_p[ii] * xres * px1_pas);
								i++;
							}
							break;
						}
						case FLOAT_DT:
						{
							float *px1_p = buf_px1;
							for (ii = x_offset_start; ii < x_offset_end; ii++) {
								SDATA(data->ew, i, j, k) = (float) (px1_p[ii] * xres * px1_pas);
								i++;
							}
							break;
						}
						default:
							error++;
							break;
						}
						j++;
					}
				} else if (config == PLANARCONFIG_SEPARATE) {
					fprintf(stderr, "px1 file is layered: case not coded yet\n");
					error++;
				}
				_TIFFfree(buf_px1);
				TIFFClose(tif_px1);
			} else {
				fprintf(stderr, "error opening px1 file\n");
				error++;
			}
		}
	}
	return error;
}

int load_px2(struct all_data *restrict data, const struct domain_info *domain, const int ind_block_x, const int ind_block_y, const float yres, const float px2_pas, const int numthreads)
{
	int error = 0;
	int num_block = ind_block_y * domain->nb_block_x + ind_block_x;
#pragma omp parallel for schedule(runtime) default(none) shared(num_block, error, data, domain, stderr) num_threads(numthreads)
	for (int k = 0; k < data->ew->nlayer; k++) {
		if (data->px2_is_tiled[k]) {
			TIFF* tif_px2 = TIFFOpen(data->px2_paths[k], "r");
			if (tif_px2 != NULL) {
			        uint32 imagewidth, imagelength;
			        uint32 tilewidth, tilelength;
				uint16 config, datatype;
				uint32 i_tile, j_tile, i_buf, j_buf;
				tdata_t buf_px2;

				int i_init, j_init;				
				TIFFGetField(tif_px2, TIFFTAG_IMAGEWIDTH, &imagewidth);
				TIFFGetField(tif_px2, TIFFTAG_IMAGELENGTH, &imagelength);
				TIFFGetField(tif_px2, TIFFTAG_TILEWIDTH, &tilewidth);
				TIFFGetField(tif_px2, TIFFTAG_TILELENGTH, &tilelength);
				TIFFGetField(tif_px2, TIFFTAG_PLANARCONFIG, &config);
				TIFFGetField(tif_px2, TIFFTAG_DATATYPE, &datatype);
				// printf("PX2 type: %d -> %s\n",datatype, datatype==INT16_DT?"int 16":datatype==FLOAT_DT?"float 32":"unknown");

				buf_px2 = _TIFFmalloc(TIFFTileSize(tif_px2));
				if (config == PLANARCONFIG_CONTIG || 1) {
					int first_x = domain->blocks[num_block].first_x;
					int first_y = domain->blocks[num_block].first_y;
					int last_x = domain->blocks[num_block].last_x;
					int last_y = domain->blocks[num_block].last_y;
					long y_offset_start = domain->blocks[num_block].y_offset_start;
					long y_offset_end = domain->blocks[num_block].y_offset_end;
					j_init = 0;
					if (!first_y) {
						y_offset_start -= data->ew->row_padding;
						j_init -= data->ew->row_padding;
					}
					if (!last_y) {
						y_offset_end += data->ew->row_padding;
					}
					long tile_start_y = y_offset_start / tilelength;
					long tile_end_y = y_offset_end / tilelength;
					long offset_start_intile_y = y_offset_start - tile_start_y * tilelength;
					long offset_end_intile_y = y_offset_end - tile_end_y * tilelength;

					long x_offset_start = domain->blocks[num_block].x_offset_start;
					long x_offset_end = domain->blocks[num_block].x_offset_end;
					i_init = 0;
					if (!first_x) {
						x_offset_start -= data->ew->col_padding;
						i_init -= data->ew->col_padding;
					}
					if (!last_x) {
						x_offset_end += data->ew->col_padding;
					}
					long tile_start_x = x_offset_start / tilewidth;
					long tile_end_x = x_offset_end / tilewidth;
					long offset_start_intile_x = x_offset_start - tile_start_x * tilewidth;
					long offset_end_intile_x = x_offset_end - tile_end_x * tilewidth;
					// assumption : image is not compressed : random access is not possible for compressed tiff
					// todo : test compression
					long shift_tile_y = 0;
					for (j_tile = tile_start_y; j_tile <= tile_end_y; j_tile++) {
						long offset_start_buf_y = (j_tile == tile_start_y) ? offset_start_intile_y : 0;
						long offset_end_buf_y = (j_tile == tile_end_y) ? offset_end_intile_y : tilelength;

						long shift_tile_x = 0;
						for (i_tile = tile_start_x; i_tile <= tile_end_x; i_tile++) {
							long offset_start_buf_x = (i_tile == tile_start_x) ? offset_start_intile_x : 0;
							long offset_end_buf_x = (i_tile == tile_end_x) ? offset_end_intile_x : tilewidth;
							TIFFReadTile(tif_px2, buf_px2, i_tile * tilewidth, j_tile * tilelength, 0, 0);
							int inbuf_j = 0;
							switch(datatype) {
							case INT16_DT:
							{
								int16_t *px2_p = buf_px2;
								for (j_buf = offset_start_buf_y; j_buf < offset_end_buf_y; j_buf++) {
									int inbuf_i = 0;
									for (i_buf = offset_start_buf_x; i_buf < offset_end_buf_x; i_buf++) {
										int i = i_init + shift_tile_x + inbuf_i;
										int j = j_init + shift_tile_y + inbuf_j;
										SDATA(data->ns, i, j, k) = (float) (px2_p[j_buf * tilewidth + i_buf] * yres * px2_pas);
										inbuf_i++;
									}
									inbuf_j++;
								}
								break;
							}
							case FLOAT_DT:
							{
								float *px2_p = buf_px2;
								for (j_buf = offset_start_buf_y; j_buf < offset_end_buf_y; j_buf++) {
									int inbuf_i = 0;
									for (i_buf = offset_start_buf_x; i_buf < offset_end_buf_x; i_buf++) {
										int i = i_init + shift_tile_x + inbuf_i;
										int j = j_init + shift_tile_y + inbuf_j;
										SDATA(data->ns, i, j, k) = (float) (px2_p[j_buf * tilewidth + i_buf] * yres * px2_pas);
										inbuf_i++;
									}
									inbuf_j++;
								}
								break;
							}
							default:
								error++;
								break;
							}
							shift_tile_x += (tilewidth - offset_start_buf_x);
						}
						shift_tile_y += (tilelength - offset_start_buf_y);
					}
				} else if (config == PLANARCONFIG_SEPARATE) {
					fprintf(stderr, "px2 file is layered: case not coded yet\n");
					error++;
				}
				_TIFFfree(buf_px2);
				TIFFClose(tif_px2);
			} else {
				fprintf(stderr, "error opening px2 file\n");
				error++;
			}
		} else {
			TIFF* tif_px2 = TIFFOpen(data->px2_paths[k], "r");
			if (tif_px2 != NULL) {
				uint32 imagelength;
				uint16 config, datatype;
				uint32 ii, jj;
				tdata_t buf_px2;

				int i, j;				
				TIFFGetField(tif_px2, TIFFTAG_IMAGELENGTH, &imagelength);
				TIFFGetField(tif_px2, TIFFTAG_PLANARCONFIG, &config);
				TIFFGetField(tif_px2, TIFFTAG_DATATYPE, &datatype);
				// printf("PX2 type: %d -> %s\n",datatype, datatype==INT16_DT?"int 16":datatype==FLOAT_DT?"float 32":"unknown");

				buf_px2 = _TIFFmalloc(TIFFScanlineSize(tif_px2));
				if (config == PLANARCONFIG_CONTIG || 1) {
					int first_x = domain->blocks[num_block].first_x;
					int first_y = domain->blocks[num_block].first_y;
					int last_x = domain->blocks[num_block].last_x;
					int last_y = domain->blocks[num_block].last_y;
					int y_offset_start = domain->blocks[num_block].y_offset_start;
					int y_offset_end = domain->blocks[num_block].y_offset_end;
					j = 0;
					if (!first_y) {
						y_offset_start -= data->ew->row_padding;
						j -= data->ew->row_padding;
					}
					if (!last_y) {
						y_offset_end += data->ew->row_padding;
					}
					assert(y_offset_end <= imagelength);
					// assumption : image is not compressed : random access is not possible for compressed tiff
					// todo : test compression
					for (jj = y_offset_start; jj < y_offset_end; jj++) {
						TIFFReadScanline(tif_px2, buf_px2, jj, 0);
						int x_offset_start = domain->blocks[num_block].x_offset_start;
						int x_offset_end = domain->blocks[num_block].x_offset_end;
						i = 0;
						if (!first_x) {
							x_offset_start -= data->ew->col_padding;
							i -= data->ew->col_padding;
						}
						if (!last_x) {
							x_offset_end += data->ew->col_padding;
						}
						switch(datatype) {
						case INT16_DT:
						{
							int16_t *px2_p = buf_px2;
							for (ii = x_offset_start; ii < x_offset_end; ii++) {
								SDATA(data->ns, i, j, k) = (float) (px2_p[ii] * yres * px2_pas);
								i++;
							}
							break;
						}
						case FLOAT_DT:
						{
							float *px2_p = buf_px2;
							for (ii = x_offset_start; ii < x_offset_end; ii++) {
								SDATA(data->ns, i, j, k) = (float) (px2_p[ii] * yres * px2_pas);
								i++;
							}
							break;
						}
						default:
							error++;
							break;
						}
						j++;
					}
				} else if (config == PLANARCONFIG_SEPARATE) {
					fprintf(stderr, "px2 file is layered: case not coded yet\n");
					error++;
				}
				_TIFFfree(buf_px2);
				TIFFClose(tif_px2);
			} else {
				fprintf(stderr, "error opening px2 file\n");
				error++;
			}
		}
	}
	return error;
}

int load_cc(struct all_data *restrict data, const struct domain_info *domain, const int ind_block_x, const int ind_block_y, const int numthreads)
{
	int error = 0;
	int num_block = ind_block_y * domain->nb_block_x + ind_block_x;
#pragma omp parallel for schedule(runtime) default(none) shared(num_block, error, data, domain, stderr) num_threads(numthreads)
	for (int k = 0; k < data->ew->nlayer; k++) {
		if (data->cc_is_tiled[k]) {
			TIFF* tif_cc = TIFFOpen(data->cc_paths[k], "r");
			if (tif_cc != NULL) {
			        uint32 imagewidth, imagelength;
			        uint32 tilewidth, tilelength;
				uint16 config;
				uint32 i_tile, j_tile, i_buf, j_buf;
				tdata_t buf_cc;

				int i_init, j_init;				
				TIFFGetField(tif_cc, TIFFTAG_IMAGEWIDTH, &imagewidth);
				TIFFGetField(tif_cc, TIFFTAG_IMAGELENGTH, &imagelength);
				TIFFGetField(tif_cc, TIFFTAG_TILEWIDTH, &tilewidth);
				TIFFGetField(tif_cc, TIFFTAG_TILELENGTH, &tilelength);
				TIFFGetField(tif_cc, TIFFTAG_PLANARCONFIG, &config);

				buf_cc = _TIFFmalloc(TIFFTileSize(tif_cc));
				if (config == PLANARCONFIG_CONTIG || 1) {
					int first_x = domain->blocks[num_block].first_x;
					int first_y = domain->blocks[num_block].first_y;
					int last_x = domain->blocks[num_block].last_x;
					int last_y = domain->blocks[num_block].last_y;
					long y_offset_start = domain->blocks[num_block].y_offset_start;
					long y_offset_end = domain->blocks[num_block].y_offset_end;
					j_init = 0;
					if (!first_y) {
						y_offset_start -= data->ew->row_padding;
						j_init -= data->ew->row_padding;
					}
					if (!last_y) {
						y_offset_end += data->ew->row_padding;
					}
					long tile_start_y = y_offset_start / tilelength;
					long tile_end_y = y_offset_end / tilelength;
					long offset_start_intile_y = y_offset_start - tile_start_y * tilelength;
					long offset_end_intile_y = y_offset_end - tile_end_y * tilelength;

					long x_offset_start = domain->blocks[num_block].x_offset_start;
					long x_offset_end = domain->blocks[num_block].x_offset_end;
					i_init = 0;
					if (!first_x) {
						x_offset_start -= data->ew->col_padding;
						i_init -= data->ew->col_padding;
					}
					if (!last_x) {
						x_offset_end += data->ew->col_padding;
					}
					long tile_start_x = x_offset_start / tilewidth;
					long tile_end_x = x_offset_end / tilewidth;
					long offset_start_intile_x = x_offset_start - tile_start_x * tilewidth;
					long offset_end_intile_x = x_offset_end - tile_end_x * tilewidth;
					// assumption : image is not compressed : random access is not possible for compressed tiff
					// todo : test compression
					long shift_tile_y = 0;
					for (j_tile = tile_start_y; j_tile <= tile_end_y; j_tile++) {
						long offset_start_buf_y = (j_tile == tile_start_y) ? offset_start_intile_y : 0;
						long offset_end_buf_y = (j_tile == tile_end_y) ? offset_end_intile_y : tilelength;

						long shift_tile_x = 0;
						for (i_tile = tile_start_x; i_tile <= tile_end_x; i_tile++) {
							long offset_start_buf_x = (i_tile == tile_start_x) ? offset_start_intile_x : 0;
							long offset_end_buf_x = (i_tile == tile_end_x) ? offset_end_intile_x : tilewidth;
							TIFFReadTile(tif_cc, buf_cc, i_tile * tilewidth, j_tile * tilelength, 0, 0);
							CC_DATA_TYPE *cc_p = buf_cc;
							int inbuf_j = 0;
							for (j_buf = offset_start_buf_y; j_buf < offset_end_buf_y; j_buf++) {
								int inbuf_i = 0;
								for (i_buf = offset_start_buf_x; i_buf < offset_end_buf_x; i_buf++) {
									int i = i_init + shift_tile_x + inbuf_i;
									int j = j_init + shift_tile_y + inbuf_j;
									SDATA(data->cc, i, j, k) = (float) (cc_p[j_buf * tilewidth + i_buf] - ZERO_CC_IN_TIF )/ (float)(255 - ZERO_CC_IN_TIF);
									inbuf_i++;
								}
								inbuf_j++;
							}
							shift_tile_x += (tilewidth - offset_start_buf_x);
						}
						shift_tile_y += (tilelength - offset_start_buf_y);
					}
				} else if (config == PLANARCONFIG_SEPARATE) {
					fprintf(stderr, "cc file is layered: case not coded yet\n");
					error++;
				}
				_TIFFfree(buf_cc);
				TIFFClose(tif_cc);
			} else {
				fprintf(stderr, "error opening cc file\n");
				error++;
			}
		} else {
			TIFF* tif_cc = TIFFOpen(data->cc_paths[k], "r");
			if (tif_cc != NULL) {
				uint32 imagelength;
				uint16 config;
				uint32 ii, jj;
				tdata_t buf_cc;

				int i, j;				
				TIFFGetField(tif_cc, TIFFTAG_IMAGELENGTH, &imagelength);
				TIFFGetField(tif_cc, TIFFTAG_PLANARCONFIG, &config);

				buf_cc = _TIFFmalloc(TIFFScanlineSize(tif_cc));
				if (config == PLANARCONFIG_CONTIG || 1) {
					int first_x = domain->blocks[num_block].first_x;
					int first_y = domain->blocks[num_block].first_y;
					int last_x = domain->blocks[num_block].last_x;
					int last_y = domain->blocks[num_block].last_y;
					int y_offset_start = domain->blocks[num_block].y_offset_start;
					int y_offset_end = domain->blocks[num_block].y_offset_end;
					j = 0;
					if (!first_y) {
						y_offset_start -= data->ew->row_padding;
						j -= data->ew->row_padding;
					}
					if (!last_y) {
						y_offset_end += data->ew->row_padding;
					}
					assert(y_offset_end <= imagelength);
					// assumption : image is not compressed : random access is not possible for compressed tiff
					// todo : test compression
					for (jj = y_offset_start; jj < y_offset_end; jj++) {
						TIFFReadScanline(tif_cc, buf_cc, jj, 0);
						CC_DATA_TYPE *cc_p = buf_cc;
						cc_p = buf_cc;
						int x_offset_start = domain->blocks[num_block].x_offset_start;
						int x_offset_end = domain->blocks[num_block].x_offset_end;
						i = 0;
						if (!first_x) {
							x_offset_start -= data->ew->col_padding;
							i -= data->ew->col_padding;
						}
						if (!last_x) {
							x_offset_end += data->ew->col_padding;
						}
						for (ii = x_offset_start; ii < x_offset_end; ii++) {
							SDATA(data->cc, i, j, k) = (float) (cc_p[ii] - ZERO_CC_IN_TIF )/ (float)(255 - ZERO_CC_IN_TIF);;
							i++;
						}
						j++;
					}
				} else if (config == PLANARCONFIG_SEPARATE) {
					fprintf(stderr, "cc file is layered: case not coded yet\n");
					error++;
				}
				_TIFFfree(buf_cc);
				TIFFClose(tif_cc);
			} else {
				fprintf(stderr, "error opening cc file\n");
				error++;
			}
		}
	}
	return error;
}

int get_files_from_cmd_line(struct all_data *data, const int nlayer, const int argc, char *const *argv, const int opt_ind)
{
	int num_layer, i, nrow, ncol;

	for (i = optind; i < argc; i++) {
		num_layer = (i-optind)/NUM_FILES;
		if ((i-optind)%NUM_FILES == 0) {
			strcpy(data->px1_paths[num_layer], argv[i]);
			get_image_info(argv[i], &nrow, &ncol, &data->px1_is_tiled[num_layer]);
		}
		if ((i-optind)%NUM_FILES == 1) {
			strcpy(data->px2_paths[num_layer], argv[i]);
			get_image_info(argv[i], &nrow, &ncol, &data->px2_is_tiled[num_layer]);
		}
		if ((i-optind)%NUM_FILES == 2) {
			strcpy(data->cc_paths[num_layer], argv[i]);
			get_image_info(argv[i], &nrow, &ncol, &data->cc_is_tiled[num_layer]);
		}
		if ((i-optind)%NUM_FILES == 3) {
			data->duration[num_layer] = atof(argv[i]);
		}
	}
	if ((i-optind) != (nlayer * NUM_FILES )) return 1;
	return 0;
}

int search_files_in_folder(struct all_data *data, const int nrow, const int ncol, const int padding, const struct parameters *param)
{
	int valret;
	pcre *re;
	const char *error;
	int erroffset;

	re = pcre_compile("(20[0-9]{6})[0-9]{7}.*(20[0-9]{6})[0-9]{7}", 0, &error, &erroffset, NULL);
	if (re == NULL) {
		printf("PCRE compilation failed at offset %d: %s\n", erroffset, error);
		return 1;
	}
	valret = search_files_impl(data, param->max_num_layer, param->infolder, "tif", param->only_multitemp, 0, -1, re);
	pcre_free(re);

	return valret;
}

int search_files_impl(struct all_data *data, const int max_num_layer, const char *dirname, const char* ext, const int only_multitemp, const int level, int num_layer, const pcre *re)
{
	DIR *d;
	struct dirent *dir;
	char childname[PATH_MAX];
	int res, nrow, ncol;

	// Assumption made : each px1, px2 and Correl files appears exactly 1 time in each 2-dated dir
	if (level > MAX_FOLDER_DEPTH) return 0;
	d = opendir(dirname);
	if (d) {
		while ((dir = readdir(d)) != NULL) {
			if ((res = is_directory(get_full_name(childname, dirname, dir->d_name))) > 0) {
				if (strcmp(dir->d_name, ".") && strcmp(dir->d_name, "..")) {
					struct my_date date1;
					struct my_date date2;
					int res;
					res = extract_dates(dir->d_name, &date1, &date2, re);
					if (res == 0) {
						if (only_multitemp && same_date(&date1, &date2)) continue;
						num_layer++;
						if (num_layer > (data->ew->nlayer) - 1) return 0;
						if (max_num_layer > 0 && num_layer > (max_num_layer - 1)) return 0;

						memcpy(&data->dates_1[num_layer], &date1, sizeof(struct my_date)); 
						memcpy(&data->dates_2[num_layer], &date2, sizeof(struct my_date)); 

						search_files_impl(data, max_num_layer, childname, ext, only_multitemp, level + 1, num_layer, re);
					} else {
						fprintf(stderr, "WARNING: no date matching for dir %s\n", dir->d_name);
					}
				}
			} else if (res == -1) {
				file_err_stop(childname);
			} else if ((res = is_reg_file(childname)) > 0) {
				if (strncmp(ext, get_filename_ext(dir->d_name), strlen(ext)) == 0) {
					char filename[PATH_MAX];
					strcpy(filename, dir->d_name);
					if (!strncmp(basename(filename), "Px1", 3)) {
						strcpy(data->px1_paths[num_layer], childname);	
						get_image_info(childname, &nrow, &ncol, &data->px1_is_tiled[num_layer]);
					} else if (!strncmp(basename(filename), "Px2", 3)) {
						strcpy(data->px2_paths[num_layer], childname);
						get_image_info(childname, &nrow, &ncol, &data->px2_is_tiled[num_layer]);
					} else if (!strncmp(basename(filename), "Correl", 6)) {
						strcpy(data->cc_paths[num_layer], childname);
						get_image_info(childname, &nrow, &ncol, &data->cc_is_tiled[num_layer]);
					}
				}
			}
		}
		closedir(d);
	} else {
		file_err_stop(dirname);
	}
	return 0;
}

int write_partial_results(struct image_f *outmat, struct all_data *restrict data, struct domain_info *domain, const int i_blck, const int j_blck, const int numthreads, const int mean_mat)
{
	int block_width = domain->blocks[j_blck * domain->nb_block_x + i_blck].size_x;
	int block_height = domain->blocks[j_blck * domain->nb_block_x + i_blck].size_y;

	long offset_start = (j_blck * outmat->ncol * domain->size_block_y + i_blck * domain->size_block_x);

#pragma omp parallel for schedule(runtime) default(none) shared(block_width, block_height, offset_start, outmat, data) num_threads(numthreads)
	for (int j = 0; j < block_height; j++) {
		long offset = offset_start + j * outmat->ncol;
		for (int i = 0; i < block_width; i++) {
			if (mean_mat == 1)
				outmat->data[offset++] = IDATA_RAW(data->mean_magn, i, j);
			else if (mean_mat == 2)
				outmat->data[offset++] = (float) IDATA(data->na_count, i, j);
			else if (mean_mat == 3)
				outmat->data[offset++] = (float) IDATA_RAW(data->mean_veloc_ew, i, j);
			else if (mean_mat == 4)
				outmat->data[offset++] = (float) IDATA_RAW(data->mean_veloc_ns, i, j);
			else
				outmat->data[offset++] = IDATA_RAW(data->out_mat, i, j);
		}
	}
	return 0;
}

int write_results_in_tiff(char* fileout, struct image_f *outmat, struct all_data *restrict data, struct domain_info *domain)
{
	TIFF* tif = TIFFOpen(fileout, "w");
	float *buf_out;
	buf_out = (float *) malloc(outmat->ncol * sizeof(float));
	tdata_t p_buf = buf_out;
	int offset = 0;

	if (tif == NULL) {
		fprintf(stderr, "Unable to open %s for writing\n", fileout);
		return 1;
	}
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, outmat->ncol);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, outmat->nrow);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(tif, TIFFTAG_DATATYPE, FLOAT_DT);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
	TIFFSetField(tif, TIFFTAG_ORIENTATION, (int)ORIENTATION_TOPLEFT);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG); // 1
	// TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_CCITTRLE);
	TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	for (int irow = 0; irow < outmat->nrow; irow++) {
		for (int i=0; i < outmat->ncol; i++) {
			buf_out[i] = (float) (outmat->data[offset++]);
		}
		TIFFWriteScanline(tif, p_buf, irow, 0);
	}
	TIFFWriteDirectory(tif);
	TIFFClose(tif);
	free(buf_out);
	return 0;
}

int get_image_info(const char* file, int* nrow, int* ncol, int* is_tiled)
{
	TIFF* tif = TIFFOpen(file, "r");
    	if (tif) {
    		uint32 image_width;
    		uint32 image_length;
    		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &image_width);
    		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &image_length);
    		*nrow = image_length;
    		*ncol = image_width;
    		*is_tiled = TIFFIsTiled(tif);
		TIFFClose(tif);
	} else {
		return 1;
	}
	return 0;
}

int create_dir (const char *name)
{
	int status;
	status = mkdir(name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (status == -1) {
		if (errno == EEXIST) {
			fprintf(stderr, "\nWARNING : folder %s already exists. Existing files might be overwritten...\n", name);
			return 0;
		}
		file_err_stop(name);
	}
	return 0;
}

int count_layers (const struct parameters *param)
{
	int valret;
	pcre *re;
	const char *error;
	int erroffset;

	re = pcre_compile("(20[0-9]{6})[0-9]{7}.*(20[0-9]{6})[0-9]{7}", 0, &error, &erroffset, NULL);
	if (re == NULL) {
		printf("PCRE compilation failed at offset %d: %s\n", erroffset, error);
		return 1;
	}
	valret = count_layers_impl(param->max_num_layer, param->infolder, param->only_multitemp, 0, re);
	pcre_free(re);

	return valret;
}

int count_layers_impl (const int max_num_layer, const char *dirname, const int only_multitemp, const int level, const pcre *re)
{
	DIR *d;
	struct dirent *dir;
	int sum = 0;
	char childname[PATH_MAX];
	int res;

	if (level > MAX_FOLDER_DEPTH) return 0;
	d = opendir(dirname);
	if (d) {
		while ((dir = readdir(d)) != NULL) {
			if ((res = is_directory(get_full_name(childname, dirname, dir->d_name))) > 0) {
				if (strcmp(dir->d_name, ".") && strcmp(dir->d_name, "..")) {
					struct my_date date1;
					struct my_date date2;
					int res;
					res = extract_dates(dir->d_name, &date1, &date2, re);
					if (res == 0) {
						if (only_multitemp && same_date(&date1, &date2)) continue;
						sum++;
						sum += count_layers_impl(max_num_layer, childname, only_multitemp, level+1, re);
					} else {
						fprintf(stderr, "WARNING: no date matching for dir %s\n", dir->d_name);
					}
				}
			} else if (res == -1) {
				file_err_stop(childname);
			}
			// DGN quick hack: carefull with this option (-k) and MAX_FOLDER_DEPTH > 1 !!! 
			if (max_num_layer > 0 && sum == max_num_layer) return sum;
		}
		closedir(d);
	} else {
		file_err_stop(dirname);
	}
	return sum;
}

// pattern matching done with libpcre (Perl-style regular expression engine)
int extract_dates (const char* dirname, struct my_date* date1, struct my_date* date2, const pcre *re)
{
	int length;
	int rc, i;
	int valret;
	char date_str[9];
	int ovector[OVECCOUNT];

	length = strlen(dirname);
	rc = pcre_exec(re, NULL, dirname, length, 0, 0, ovector, OVECCOUNT);
	if (rc < 0) {
		switch (rc) {
	    	case PCRE_ERROR_NOMATCH:
	    		// No match
	    		valret = 1;
		default:
			// Matching error
			valret = 2;
		}
	} else {
		// Match succeeded
		for (i = 1; i < rc; i++) { // trick : start to 1 to avoid the complete matching : we just want the 2 substrings matched, not the first all one
			const char *substring_start = dirname + ovector[2*i];
			int substring_length = ovector[2*i+1] - ovector[2*i];

			strncpy(date_str, substring_start, substring_length);
			date_str[8] = '\0';
			switch(i) {
			case 1:
				get_date(date1, date_str);
				break;
			case 2:
				get_date(date2, date_str);
				break;
			}
		}	
		valret = 0;
	}
	return valret;
}

const char *get_filename_ext (const char *filename)
{
	const char *dot = strrchr(filename, '.');
	if(!dot || dot == filename) return "";
	return dot + 1;
}

int is_directory(const char* path)
{
	struct stat buf;
	int res;

	res = stat(path, &buf);
	if (res == 0) {
		return S_ISDIR(buf.st_mode)?1:0;
	} else {
		return -1;
	}
}

int is_reg_file(const char* path)
{
	struct stat buf;

	int res;

	res = stat(path, &buf);
	if (res == 0) {
		return S_ISREG(buf.st_mode)?1:0;
	} else {
		return -1;
	}
}

char *get_full_name(char *fullname, const char *parent, const char *child)
{
	strcpy(fullname, parent);
	strcat(fullname, "/");
	strcat(fullname, child);
	return fullname;
}

void file_err_stop(const char* filename)
{
	char err_str[PATH_MAX + 20];
	sprintf(err_str, "Error accessing %s", filename);
	perror(err_str);
	abort();
}

int write_partial_results_in_bin_file(FILE *fileout, struct all_data *restrict data, struct domain_info *domain, const int i, const int j, const int nrow, const int ncol, const int mean_mat)
{
	int block_width = domain->blocks[j * domain->nb_block_x + i].size_x;
	int block_height = domain->blocks[j * domain->nb_block_x + i].size_y;

	long offset_start = (j * ncol * domain->size_block_y + i * domain->size_block_x) * sizeof(float);
	long offset_next_line = (ncol - block_width) * sizeof(float);

	if (fseek(fileout, offset_start, SEEK_SET) < 0) return 1;
	for (int irow = 0; irow < block_height; irow++) {
		for (int icol = 0; icol < block_width; icol++) {
			if (mean_mat == 1)
				fwrite(&IDATA_RAW(data->mean_magn, icol, irow), sizeof(float), 1, fileout);		
			else if (mean_mat == 2)
				fwrite(&IDATA(data->na_count, icol, irow), sizeof(int), 1, fileout);		
			else if (mean_mat == 3)
				fwrite(&IDATA_RAW(data->mean_veloc_ew, icol, irow), sizeof(float), 1, fileout);		
			else if (mean_mat == 4)
				fwrite(&IDATA_RAW(data->mean_veloc_ns, icol, irow), sizeof(float), 1, fileout);		
			else
				fwrite(&IDATA_RAW(data->out_mat, icol, irow), sizeof(float), 1, fileout);		
		}
		if (fseek(fileout, offset_next_line, SEEK_CUR) < 0) return 1;
	}
	return 0;
}

int write_results_in_tiff_from_binfile(char* fileout, struct all_data *restrict data, struct domain_info *domain, const int ncol, const int nrow, const int mean_mat)
{
	FILE* fin;

	if (mean_mat == 1) 
		fin = fopen(MEAN_TMP_FILE, "r" );	
	else if (mean_mat == 2)
		fin = fopen(NACOUNT_TMP_FILE, "r" );
	else if (mean_mat == 3)
		fin = fopen(VELOC_EW_TMP_FILE, "r" );
	else if (mean_mat == 4)
		fin = fopen(VELOC_NS_TMP_FILE, "r" );
	else
		fin = fopen(BIN_TMP_FILE, "r" );

	TIFF* tif = TIFFOpen(fileout, "w");
	float *buf_in;
	buf_in = (float *) malloc(ncol * sizeof(float));
	
	tdata_t p_buf;
	p_buf = buf_in;

	if (tif == NULL) {
		fprintf(stderr, "Unable to open %s for writing\n", fileout);
		return 1;
	}
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, ncol);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, nrow);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(tif, TIFFTAG_DATATYPE, FLOAT_DT);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
	TIFFSetField(tif, TIFFTAG_ORIENTATION, (int)ORIENTATION_TOPLEFT);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG); // 1
	TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	for (int irow = 0; irow < nrow; irow++) {
		assert(fread(buf_in, sizeof(float), ncol, fin) == ncol);
		TIFFWriteScanline(tif, p_buf, irow, 0);
	}
	TIFFWriteDirectory(tif);
	TIFFClose(tif);
	fclose(fin);
	free(buf_in);
	return 0;
}

#ifdef FAST_HYPOT
#ifndef EVEN_FASTER
static inline float fast_sqrtf(float number)
{
// avoid gcc warnings for this well known hack
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y  = number;
	i  = * ( long * ) &y;
	i  = 0x5f375a86 - ( i >> 1 );               // Lomont's magic number
	y  = * ( float * ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
#ifdef SECOND_NEWTON_PASS
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration (can be removed)
#endif
#pragma GCC diagnostic pop
	return 1/y;
}

float fast_hypotf(float x, float y)
{
        /*      Computes sqrt(x*x+y*y), avoiding overflow */ 
        if (x < 0) x = -x;
        if (y < 0) y = -y;
        if (x > y) {
                float t = y;
                y = x;
                x = t;
        }
        /* sqrt(x*x+y*y) = sqrt(y*y*(x*x/(y*y)+1.0)) = y*sqrt(x*x/(y*y)+1.0) */
        if (y == 0.0) return 0.0;
        x /= y;
        return y*fast_sqrtf(x*x+1.0);
}
#else
float fast_hypotf(float x, float y)
{	// implementation of alpha max plus beta min algorithm
	// largest error: 3.96%, mean error: 2.41% 
        float min, max;
        const float alpha = 0.96043387F;
        const float beta = 0.39782473F;

        if (x < 0) x = -x;
        if (y < 0) y = -y;
        if (x > y) {
                max = x;
                min = y;
        } else {
                min = x;
                max = y;
	}
	return alpha * max + beta * min;
}
#endif
#endif

#include "func.h"
#include "misc.h"

int main (const int argc, char *const *argv)
{
	struct parameters params = {	.infolder = NULL,
					.xres = 0.5, .yres = 0.5,
					.px1_pas = 0.1, .px2_pas = 0.1,
					.ref_image = NULL,
					.neighborhood = 11,
					.min_displ = 0.15,
					.min_cc = 0.333,
					.only_multitemp = 1,
					.min_matches = 0.4,
					.available_GiB_mem = 2.f,
					.output_dir = NULL,
					.output_mean = 0,
					.output_nacount = 0,
					.output_veloc = 0,
					.max_num_layer = 0,
					.nlayer = 0 };

	struct all_data *data[2];
	struct image_f *outmat_results;
	struct image_f *outmat_mean;
	struct image_f *outmat_nacount;
	struct image_f *outmat_veloc_ew;
	struct image_f *outmat_veloc_ns;
	struct domain_info *domain;

	struct all_time timings;
	int nrow, ncol, padding, is_tiled, path_opt_ind;
	int only_one_block = 0;
	int output_in_file;

  	char outdir[2*PATH_MAX];
  	char out_pattern[2*PATH_MAX];
  	char out_tif[2*PATH_MAX];
  	char out_mean[2*PATH_MAX];
  	char out_nacount[2*PATH_MAX];
  	char out_veloc_ew[2*PATH_MAX];
  	char out_veloc_ns[2*PATH_MAX];
  	FILE *binfile;
  	FILE *binfile_mean;
  	FILE *binfile_nacount;
  	FILE *binfile_veloc_ew;
  	FILE *binfile_veloc_ns;

  	int nb_threads, nthreads_read, nthreads_comp;
	int comp_set = 1;
	int write_set = 0;

	reset_timer(&timings);
	CLOCK_START(timings.glob);

	if (VERBOSE) print_switches();
	
	/* get options from command line to fill params */
	if (get_cmd_line_options(&params, argc, argv, &path_opt_ind)) {
		print_usage(argv);
		return 1;
	}

	padding = (params.neighborhood - 1)/2;
	if (params.nlayer == 0) {
		/* count the number of layer in the underlying folders of params.infolder */
		params.nlayer = count_layers(&params);
		if (VERBOSE) {
			fprintf(stdout, "***************** LAYERS ************************\n");
			printf("%d layers found in %s\n", params.nlayer, params.infolder);
			printf("*************************************************\n\n");
		}
	}
	/* disable libtiff warning handler so that it stops complaining about the non standard format of our input files */
  	TIFFSetWarningHandler(0);

  	/* get some info about the image stack from the reference image: nrow, ncol, and if internal layout is tiled or striped */
  	if (get_image_info(params.ref_image, &nrow, &ncol, &is_tiled)) {
  		fprintf(stderr, "failed to open file %s\n", params.ref_image);
  		perror("reason");
		return 1;
  	}


	if (VERBOSE > 1) {
		fprintf(stdout, "****************** IMAGE ************************\n");
		printf("%i x %i\n", ncol, nrow);
		printf("%s\n", is_tiled?"TILED":"STRIPPED");
		printf("*************************************************\n\n");
	}

	/* creating dir & names for outputs */
	if (params.output_dir == NULL) {
		strcpy(outdir, ".");
	} else {
		strcpy(outdir, params.output_dir);
		if (mkdir_rec(outdir) != 0) {
			fprintf(stderr, "Check your output dir: %s\naborting.\n", outdir);
			_exit(1);
		}
	}
	strcpy(out_pattern, outdir);
	strcat(out_pattern, "/outmat");
	
	strcpy(out_tif, out_pattern);
	strcat(out_tif, "_vc.tif");
	
	strcpy(out_mean, out_pattern);
	strcat(out_mean, "_mean.tif");
	
	strcpy(out_nacount, out_pattern);
	strcat(out_nacount, "_nacount.tif");
	
	strcpy(out_veloc_ew, out_pattern);
	strcat(out_veloc_ew, "_veloc_ew.tif");
	
	strcpy(out_veloc_ns, out_pattern);
	strcat(out_veloc_ns, "_veloc_ns.tif");

#ifdef _OPENMP
  	/* Set 2 teams of threads : 1 set for loading data, the other set for computing*/
  	omp_set_nested(1);
  	/* check if OMP_NUM_THREADS env var is defined */
  	char * omp_num_threads = getenv("OMP_NUM_THREADS");
  	if (omp_num_threads == NULL) {
  		// take all the threads available (maximum of node capacity)
  		nb_threads = omp_get_num_procs(); 
  	} else {
  		// take the number of thread defined by OMP_NUM_THREADS
  		nb_threads = atoi(omp_num_threads);
  	}
  	// temporary cooking ... need to be tuned according to the machine I/O vs computation performance
/*	nthreads_read = 6;
	nthreads_comp = nb_threads - nthreads_read;*/
  	if (nb_threads >= params.nlayer) {
		nthreads_read = (params.nlayer / 2)>=1?(params.nlayer / 2):1;
		nthreads_comp = nb_threads - nthreads_read;
	} else if (nb_threads >= 4) {
		nthreads_read = (int) ceilf(nb_threads / 4.0);
		// nthreads_read = (int) ceilf(nb_threads / 2.0);
		nthreads_comp = nb_threads - nthreads_read;
	} else {
		nthreads_read = 1;
		nthreads_comp = nb_threads - nthreads_read;
	}
#else
	nb_threads = 1;
	nthreads_read = 1;
	nthreads_comp = 1;
#endif
	/* if the problem doesn't fit in memory, cut the domain into blocks */
	if (VERBOSE) fprintf(stdout, "**************** AUTO CONF **********************\n");

	if (init_domain_decomposition(&params, &domain, params.available_GiB_mem, ncol, nrow, params.nlayer, padding, padding, &output_in_file)) {
		fprintf(stderr, "init_domain_decomposition failed !\n");
		return 1;
	}
	if (output_in_file) {
		binfile = fopen(BIN_TMP_FILE, "w");
		if (binfile == NULL) file_err_stop(BIN_TMP_FILE);
		if (params.output_mean) {
			binfile_mean = fopen(MEAN_TMP_FILE, "w");
			if (binfile_mean == NULL) file_err_stop(MEAN_TMP_FILE);
		}
		if (params.output_nacount) {
			binfile_nacount = fopen(NACOUNT_TMP_FILE, "w");
			if (binfile_nacount == NULL) file_err_stop(NACOUNT_TMP_FILE);
		}
		if (params.output_veloc) {
			binfile_veloc_ew = fopen(VELOC_EW_TMP_FILE, "w");
			if (binfile_veloc_ew == NULL) file_err_stop(VELOC_EW_TMP_FILE);
			binfile_veloc_ns = fopen(VELOC_NS_TMP_FILE, "w");
			if (binfile_veloc_ns == NULL) file_err_stop(VELOC_NS_TMP_FILE);
		}
	}
	if (domain->nb_block_x == 1 && domain->nb_block_y == 1) {
		only_one_block = 1;
		nthreads_read = nb_threads;
		nthreads_comp = nb_threads;
	}
	if (VERBOSE) {
		// fprintf(stdout, "**************** AUTO CONF **********************\n");
		printf("domain decomposition:\r\t\t\t%d x %d blocks\n", domain->nb_block_x, domain->nb_block_y);
		printf("temporary results:\r\t\t\t%s\n", output_in_file?"on disk":"in memory");
		printf("*************************************************\n\n");
	}
	if (VERBOSE > 2) printf("******************* OpenMP **********************\n");
#ifdef _OPENMP
	if (VERBOSE > 2) {
		printf("nb threads:\r\t\t\t%d\n", nb_threads);
		printf("\tfor loading:\r\t\t\t%d\n", nthreads_read);
		printf("\tfor computing:\r\t\t\t%d\n", nthreads_comp);
	}
#else
	if (VERBOSE > 2) printf("Compiled without OpenMP\n");
#endif
	if (VERBOSE > 2) printf("*************************************************\n\n");

	printf("running ...\n");

	/* allocate memory for data structure */
	data[0] = allocate_all_data(domain->size_block_y, domain->size_block_x, params.nlayer, padding, padding);
	if (!only_one_block) data[1] = allocate_all_data(domain->size_block_y, domain->size_block_x, params.nlayer, padding, padding);
	if (!output_in_file) {
		outmat_results = allocate_image_f(nrow, ncol, padding, padding);
		if (params.output_mean)
			outmat_mean = allocate_image_f(nrow, ncol, padding, padding);
		if (params.output_nacount)
			outmat_nacount = allocate_image_f(nrow, ncol, padding, padding);
		if (params.output_veloc) {
			outmat_veloc_ew = allocate_image_f(nrow, ncol, padding, padding);
			outmat_veloc_ns = allocate_image_f(nrow, ncol, padding, padding);
		}
	}

	/* search/get file names and store them in data */
	if (params.infolder != NULL) {
		if (search_files_in_folder(data[0], domain->size_block_y, domain->size_block_x, padding, &params)) {
			fprintf(stderr, "search_files failed!\n");
			return 1;
		}
		if (!only_one_block) {
			if (search_files_in_folder(data[1], domain->size_block_y, domain->size_block_x, padding, &params)) {
				fprintf(stderr, "search_files failed!\n");
				return 1;
			}
		}
	} else {
		if (get_files_from_cmd_line(data[0], params.nlayer, argc, argv, path_opt_ind)) {
			fprintf(stderr, "get_files failed !\n");
			return 1;
		}
		if (!only_one_block) {
			if (get_files_from_cmd_line(data[1], params.nlayer, argc, argv, path_opt_ind)) {
				fprintf(stderr, "get_files failed !\n");
				return 1;
			}
		}
	}

	/* set arrays to zero (first touch)*/
	if (VERBOSE > 2) printf("First touch\n");
	clear_stacks_first_touch(data[0], nthreads_comp);
	if (!only_one_block) clear_stacks_first_touch(data[1], nthreads_comp);

	/* Initial load of data */
	if (VERBOSE > 2) printf("Initial load\n");
	load_stacks(data[write_set], domain, 0, 0, &params, nb_threads);


	/* Loop on blocks: out of core algorithm with read/compute overlap */
	CLOCK_START(timings.loop);
#pragma omp parallel default(none) shared(data, domain, params, is_tiled, nrow, ncol, nthreads_read, nthreads_comp, timings, write_set, comp_set, outmat_results, outmat_mean, outmat_nacount, outmat_veloc_ew, outmat_veloc_ns, output_in_file, binfile, binfile_mean, binfile_nacount, binfile_veloc_ew, binfile_veloc_ns) num_threads(2)
{
#pragma omp single
{
	/* loop along blocks */
	for (int j = 0; j < domain->nb_block_y; j++) {
		for (int i = 0; i < domain->nb_block_x; i++) {
			if (VERBOSE) printf("computing block %d / %d\n", i + j * domain->nb_block_x + 1, domain->nb_block_x * domain->nb_block_y);
			// switch data sets
			write_set = comp_set?1:0;
			comp_set = write_set?0:1;

			#pragma omp task // task preload
			{
				/* preload the next block of data */
				int inext = i;
				int jnext = j;
			CLOCK_START(timings.task_preload);
				if (++inext == domain->nb_block_x) {
					inext = 0;
					jnext++;
				}
				if (jnext < domain->nb_block_y) {
#ifdef CLEAR_IN_TASK_PRELOAD
					if (!(i==0 && j==0)) {
						CLOCK_START(timings.clear);
							clear_stacks_padding_only(data[write_set], nthreads_read);
						CLOCK_STOP(timings.clear);
					}
#endif
					CLOCK_START(timings.load);
						load_stacks(data[write_set], domain, inext, jnext, &params, nthreads_read);
					CLOCK_STOP(timings.load);
				}
			CLOCK_STOP(timings.task_preload);
			CLOCK_START(timings.preload_wait);
			}

			#pragma omp task // task process
			{
			CLOCK_START(timings.task_process);
				/* compute mean_magnitude along the time (along stack depth) */
				CLOCK_START(timings.mean);
					compute_mean_magnitude(data[comp_set], nthreads_comp);
					compute_mean_velocity(data[comp_set], nthreads_comp);
				CLOCK_STOP(timings.mean);


				/* compute persistent deformation on stack: hot spot of the code*/
				CLOCK_START(timings.process);
					process_stacks(data[comp_set], &params, nthreads_comp);
				CLOCK_STOP(timings.process);

				/* write block results in temporary binary file in ./outtmp/ */
				CLOCK_START(timings.write);
					if (output_in_file) {
						write_partial_results_in_bin_file(binfile, data[comp_set], domain, i, j, nrow, ncol, 0);
						if (params.output_mean)
							write_partial_results_in_bin_file(binfile_mean, data[comp_set], domain, i, j, nrow, ncol, 1);
						if (params.output_nacount)
							write_partial_results_in_bin_file(binfile_nacount, data[comp_set], domain, i, j, nrow, ncol, 2);
						if (params.output_veloc) {
							write_partial_results_in_bin_file(binfile_veloc_ew, data[comp_set], domain, i, j, nrow, ncol, 3);
							write_partial_results_in_bin_file(binfile_veloc_ns, data[comp_set], domain, i, j, nrow, ncol, 4);
						}
					} else {
						write_partial_results(outmat_results, data[comp_set], domain, i, j, nthreads_comp, 0);
						if (params.output_mean)
							write_partial_results(outmat_mean, data[comp_set], domain, i, j, nthreads_comp, 1);
						if (params.output_nacount)
							write_partial_results(outmat_nacount, data[comp_set], domain, i, j, nthreads_comp, 2);
						if (params.output_veloc) {
							write_partial_results(outmat_veloc_ew, data[comp_set], domain, i, j, nthreads_comp, 3);
							write_partial_results(outmat_veloc_ns, data[comp_set], domain, i, j, nthreads_comp, 4);
						}
					}
				CLOCK_STOP(timings.write);
#ifndef CLEAR_IN_TASK_PRELOAD
				if (!(i == (domain->nb_block_x - 1) && j == (domain->nb_block_y - 1))) {
					CLOCK_START(timings.clear);
					clear_stacks_padding_only(data[comp_set], nthreads_comp);
					CLOCK_STOP(timings.clear);
				}
#endif
			CLOCK_STOP(timings.task_process);
			CLOCK_START(timings.process_wait);
			}

			// synchronization
			#pragma omp taskwait
			CLOCK_STOP(timings.loop);
			CLOCK_STOP(timings.process_wait);
			CLOCK_STOP(timings.preload_wait);
			if (VERBOSE > 2) {
				printf("\tpreload:\r\t\t\t%7.2f sec %s\n", CLOCK_GET(timings.task_preload), (CLOCK_GET(timings.task_preload)>CLOCK_GET(timings.task_process))?"*":"");
				printf("\tprocess:\r\t\t\t%7.2f sec %s\n", CLOCK_GET(timings.task_process), (CLOCK_GET(timings.task_preload)>CLOCK_GET(timings.task_process))?"":"*");
				printf("\tloop:\r\t\t\t%7.2f sec\n", CLOCK_GET(timings.loop));
			}
			CLOCK_START(timings.loop);
		}
	}
}
}
	CLOCK_STOP(timings.glob);

	/* write back result to geotiff file */
	if (VERBOSE) printf("write results in tiff format\n");
	if (output_in_file) {
		write_results_in_tiff_from_binfile(out_tif, data[comp_set], domain, ncol, nrow, 0);
		fclose(binfile);
		remove(BIN_TMP_FILE);
		if (params.output_mean) {
			write_results_in_tiff_from_binfile(out_mean, data[comp_set], domain, ncol, nrow, 1);
			fclose(binfile_mean);
			remove(MEAN_TMP_FILE);
		}
		if (params.output_nacount) {
			write_results_in_tiff_from_binfile(out_nacount, data[comp_set], domain, ncol, nrow, 2);
			fclose(binfile_nacount);
			remove(NACOUNT_TMP_FILE);
		}
		if (params.output_veloc) {
			write_results_in_tiff_from_binfile(out_veloc_ew, data[comp_set], domain, ncol, nrow, 3);
			fclose(binfile_veloc_ew);
			remove(VELOC_EW_TMP_FILE);
			write_results_in_tiff_from_binfile(out_veloc_ns, data[comp_set], domain, ncol, nrow, 4);
			fclose(binfile_veloc_ns);
			remove(VELOC_NS_TMP_FILE);
		}
	} else {
		write_results_in_tiff(out_tif, outmat_results, data[comp_set], domain);
		if (params.output_mean)
			write_results_in_tiff(out_mean, outmat_mean, data[comp_set], domain);
		if (params.output_nacount)
			write_results_in_tiff(out_nacount, outmat_nacount, data[comp_set], domain);
		if (params.output_veloc) {
			write_results_in_tiff(out_veloc_ew, outmat_veloc_ew, data[comp_set], domain);
			write_results_in_tiff(out_veloc_ns, outmat_veloc_ns, data[comp_set], domain);
		}
	}

	if (VERBOSE) printf("copy metadata from ref image\n");
	if (copy_geotiff_metadata(out_tif, params.ref_image)) {
		fprintf(stderr, "Error in copy_geotiff_metadata\n");
	}
	if (params.output_mean && copy_geotiff_metadata(out_mean, params.ref_image)) {
		fprintf(stderr, "Error in copy_geotiff_metadata for mean output\n");
	}
	if (params.output_veloc && copy_geotiff_metadata(out_veloc_ew, params.ref_image)) {
		fprintf(stderr, "Error in copy_geotiff_metadata for ew veloc output\n");
	}
	if (params.output_veloc && copy_geotiff_metadata(out_veloc_ns, params.ref_image)) {
		fprintf(stderr, "Error in copy_geotiff_metadata for ns veloc output\n");
	}


	if (VERBOSE > 1) {
		print_time(&timings);
	}

	printf("done !\n");

	/* free all data structures */
	free_domain(domain);
	free_all_data(data[0]);
	if (!only_one_block) free_all_data(data[1]);
	if (!output_in_file) {
		free_image_f(outmat_results);
		if (params.output_mean)
			free_image_f(outmat_mean);
	}
	return 0;
}

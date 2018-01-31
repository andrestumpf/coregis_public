#include "misc.h"

void print_switches(void)
{
	printf("******************* SWITCHES ********************\n");
#ifdef FASTER_NACOUNT
	printf("FASTER_NACOUNT:\r\t\t\tenabled\n");
#endif
#ifdef WINDOW_SIZE
	printf("WINDOW_SIZE:\r\t\t\tenabled\r\t\t\t\t(!)\n");
	printf("\t-> WARNING !\r\t\t\tfor performance reason, use a fixed window size (%d) instead of the one passed by option -n\n\t\t\tif needed, comment the flag or modify the value in data.h and recompile\n", WINDOW_SIZE);
#endif
#ifdef FAST_HYPOT
#ifdef EVEN_FASTER
	printf("FAST_HYPOT:\r\t\t\tenabled\r\t\t\t\t(!)\n");
	printf("EVEN_FASTER:\r\t\t\tenabled\r\t\t\t\t(!)\n");
	printf("\t-> WARNING !\r\t\t\tuse a very fast sqrt implementation:\n\t\t\talpha max plus beta min algorithm:\n\t\t\t\tlargest error: 3.96%%\n\t\t\t\tmean error: 2.41%%\n");
#else
	printf("FAST_HYPOT:\r\t\t\tenabled\r\t\t\t\t(!)\n");
	printf("\t-> WARNING !\r\t\t\tuse a fast sqrt implementation:\n\t\t\tbased on Fast InvSqrt() with Chris Lomont's magic number:\n\t\t\t\tmean error: 0.175%%\n");
#ifdef SECOND_NEWTON_PASS
	printf("SECOND_NEWTON_PASS:\r\t\t\tenabled\n");
	printf("\t\t\tapply a second pass of Newton's method to increase precision of FAST_HYPOT\n");
#endif
#endif
#endif
	printf("*************************************************\n");
	printf("\n");
}

#define FASTER_NACOUNT    		/* compute na_count_blck faster (should stay on, unless you want to modify this part of code ...)			     */
// #define FAST_HYPOT			/* use fast sqrt approximation inside hypot funct -> see in function code how to increase precision    			     */
#define EVEN_FASTER


void print_time(struct all_time *timings)
{
	printf("\n");
	printf("*************** TIME MEASUREMENTS ***************\n");
#ifdef CLEAR_IN_TASK_PRELOAD
	printf("preload:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->task_preload.sum, (timings->task_preload.sum / timings->loop.sum) * 100.);
	printf("\tloads:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->load.sum, (timings->load.sum / timings->loop.sum) * 100.);
	printf("\tclear:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->clear.sum, (timings->clear.sum / timings->loop.sum) * 100.);
	printf("\twait:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->preload_wait.sum, (timings->preload_wait.sum / timings->loop.sum) * 100.);
	printf("process:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->task_process.sum, (timings->task_process.sum / timings->loop.sum) * 100.);
	printf("\tmean:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->mean.sum, (timings->mean.sum / timings->loop.sum) * 100.);
	printf("\tprocess:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->process.sum, (timings->process.sum / timings->loop.sum) * 100.);
	printf("\twrite:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->write.sum, (timings->write.sum / timings->loop.sum) * 100.);
	printf("\twait:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->process_wait.sum, (timings->process_wait.sum / timings->loop.sum) * 100.);
#else
	printf("preload:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->task_preload.sum, (timings->task_preload.sum / timings->loop.sum) * 100.);
	printf("\tloads:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->load.sum, (timings->load.sum / timings->loop.sum) * 100.);
	printf("\twait:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->preload_wait.sum, (timings->preload_wait.sum / timings->loop.sum) * 100.);
	printf("process:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->task_process.sum, (timings->task_process.sum / timings->loop.sum) * 100.);
	printf("\tmean:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->mean.sum, (timings->mean.sum / timings->loop.sum) * 100.);
	printf("\tprocess:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->process.sum, (timings->process.sum / timings->loop.sum) * 100.);
	printf("\twrite:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->write.sum, (timings->write.sum / timings->loop.sum) * 100.);
	printf("\tclear:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->clear.sum, (timings->clear.sum / timings->loop.sum) * 100.);
	printf("\twait:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->process_wait.sum, (timings->process_wait.sum / timings->loop.sum) * 100.);
#endif
	printf("-------------------------------------------------\n");
	printf("loop:\r\t\t\t%7.2f sec\r\t\t\t\t-> %4.1f%%\n", timings->loop.sum, 100.0);
	printf("global:\r\t\t\t%7.2f sec\r\t\t\t\t-> %d min %d sec\n", timings->glob.sum, (int)floor(timings->glob.sum / 60.0), (int)(timings->glob.sum - 60 * floor(timings->glob.sum / 60.0)));
	printf("*************************************************\n");
	printf("\n");
}

void reset_timer(struct all_time *timings)
{
	timings->load.sum = 0.;
	timings->mean.sum = 0.;
	timings->process.sum = 0.;
	timings->write.sum = 0.;
	timings->clear.sum = 0.;
	timings->glob.sum = 0.;
	timings->loop.sum = 0.;
	timings->task_preload.sum = 0.;
	timings->task_process.sum = 0.;
	timings->preload_wait.sum = 0.;
	timings->process_wait.sum = 0.;
}

void print_usage(char *const *argv)
{
	printf("Usage:  %s	{-x XRES} {-y YRES}\n", argv[0]);
 	printf("			{-u PX1PAS} {-v PX2PAS}\n");
	printf("			{-r REFIMAGE}\n");
	printf("			{-n NEIGHBORHOOD}\n");
	printf("			{-d MINDISPLACEMENT} {-c MINCC}\n");
 	printf("			[-m ONLYMULTITEMP]\n");
  	printf("			{-i <INFOLDER> | -l <nb_layer> <Px1_lay1_path> <Px2_lay1_path> <Correl_lay1_path> <Px1_lay2_path> ... <Correl_laynb_layer_path>}\n");
  	printf("use %s -h for detailed Information\n\n", argv[0]);
  	_exit(1);
}

void print_help(char *const *argv)
{
	printf("Usage:  %s	{-x XRES} {-y YRES}\n", argv[0]);
 	printf("			{-u PX1PAS} {-v PX2PAS}\n");
	printf("			{-r REFIMAGE}\n");
	printf("			{-n NEIGHBORHOOD}\n");
	printf("			{-d MINDISPLACEMENT} {-c MINCC}\n");
 	printf("			[-m ONLYMULTITEMP]\n");
  	printf("			{-i <INFOLDER> | -l <nb_layer> <Px1_lay1_path> <Px2_lay1_path> <Correl_lay1_path> <Px1_lay2_path> ... <Correl_laynb_layer_path>}\n");
	printf("\n");
	printf(" This programm will process a time series of correlograms resulting from MicMac using an out-of-core multi-threaded algorithm\n");
	printf("\n");
	printf(" Conventions for -i option:\n");
	printf("   - The input folder holds several subfolders of which each contains the EW & NS components + the correlation coefficient\n");
	printf("   - The input folder names should contain a string of the form YYYYMMDD from which the date of the image will be extracted automatically\n");
	printf("   - The output files are stored in a subfolder named outtmp.\n");
	printf("\n");
	printf(" Input arguments:\n");
	printf("	-h	Print this help\n");
	printf("	-e	output mean displacment magnitude in meters\n");
	printf("	-V	output mean velocity (NS + EW) in meters/day\n");
	printf("	-E	output na count\n");
	printf("	-i	Input folder holding the correlograms\n");
	printf("	-l	set the number of layer and pass all the files to be computed on the command line: after options, pass in order: paths_to: px1, px2, cc for layer 1 to n. See usage\n");
	printf("	-g	give the lower bound of available node memory in GiB to fit in memory : it is actually far better to let 1 GiB free for the system\n");
	printf("	-x	Resolution in x of the input image\n");
	printf("	-y	Resolution in y of the input image\n");
	printf("	-u	Subpixel resolution in x direction\n");
	printf("	-v	Subpixel resolution in y direction\n");
	printf("	-k	Consider only the k first layers\n");
	printf("        -r      One of the orthoimages from which the displacement field has been derrived (this serves as a reference)\n");
	printf("	-n	Neighborhood in which the vector coherence will be computed\n");
	printf("	-d	Minimum mean displacement threshold. Pixel below this value are not processed an appear as NA in the\n");
	printf("		result. Increasing the threshold can significantly speed up computation but may also comprise the risk\n");
	printf("		of not detecting real displacement below the threshold.\n");
	printf("	-c	Minimum correlation coefficient in [0;1]. Pixel in the stack with a correlation coefficient below this value\n");
	printf("		are set to NA and are consequently not considered in the multi-temporal analysis.\n");
	printf("	-o	output dir. If not specified, output in working dir\n");
	printf("	-t	[0,1], Stacking only multi-temporal correlation results if set to 1. Default is 0 where all images\n");
	printf("		will be stacked.\n");
	printf("	-m	[0...1], Fraction of matches that should exceed the minimum correlation threshold (-c). Where there are less than\n");
	printf("		m*timeSteps the output will be set to NA.\n");
	printf("\n");
	printf("Example: %s -i '/mnt/hdfs/UBO/PLEIADES/CORRELOGRAMS_2600' -x 0.5 -y 0.5 -u 0.2 -v 0.2 -r '/mnt/hdfs/UBO/PLEIADES/ORTHOIMAGES/REF_SUBSET_2600/orthoimg_phr1b_p_201406201038211_sen_948132101-002_r1c1_cropped_sub.tif' -n 5 -d 0.15 -c 0.7 -t 0 -m 0.4 -g 3\n", argv[0]);
	printf("\n");
  	_exit(1);
}

int get_cmd_line_options(struct parameters *params, const int argc, char *const *argv, int *path_opt_ind)
{
	int c;
	int definfolder = 0;
	int deflayer = 0;
	int defmultitemp = 0;
	int def_avail_mem = 0;
	int limit_num_layer = 0;
	opterr = 0;

	while ((c = getopt (argc, argv, ":g:l:i:x:y:u:v:r:n:d:k:c:t:m:heEVo:")) != -1) {
		switch (c) {
		case 'g':
			params->available_GiB_mem = atof(optarg);
			def_avail_mem = 1;
			break;
		case 'l':
			params->nlayer = atoi(optarg);
			deflayer = 1;
			break;
		case 'i':
			params->infolder = optarg;
			definfolder = 1;
			break;
		case 'x':
			params->xres = atof(optarg);
			break;
		case 'y':
			params->yres = atof(optarg);
			break;
		case 'u':
			params->px1_pas = atof(optarg);
			break;
		case 'v':
			params->px2_pas = atof(optarg);
			break;
		case 'r':
			params->ref_image = optarg;
			break;
		case 'k':
			params->max_num_layer = atoi(optarg);
			limit_num_layer = 1;
			break;
		case 'o':
			params->output_dir = optarg;
			break;
		case 'n':
			params->neighborhood = atoi(optarg);
			if (params->neighborhood < 0 || (params->neighborhood%2) == 0) {
				fprintf (stderr, "Option -n must provide an unsigned odd int\n");
				return 1;
			}
			break;
		case 'd':
			params->min_displ = atof(optarg);
			break;
		case 'c':
			params->min_cc = atof(optarg);
			break;
		case 't':
			params->only_multitemp = atoi(optarg);
			if (params->only_multitemp != 0 && params->only_multitemp != 1) {
				fprintf (stderr, "Option -t must provide an integer in [0, 1]\n");
				return 1;
			}
			defmultitemp = 1;
			break;
		case 'm':
			params->min_matches = atof(optarg);
			if (params->min_matches < 0.f || params->min_matches > 1.f) {
				fprintf (stderr, "Option -m must provide a float in [0, 1]\n");
				return 1;
			}
			break;
		case 'h':
			print_help(argv);
			break;
		case 'e':
			params->output_mean = 1;
			break;
		case 'E':
			params->output_nacount = 1;
			break;
		case 'V':
			params->output_veloc = 1;
			break;
		case ':': // missing option argument
			fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			return 1;
			break;
		case '?': // unknown option
			if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
			return 1;
		default:
			fprintf (stderr, "Problem parsing cmd line, abort\n");
			_exit(1);
		}
	}

	// check parameters coherence
	if (params->ref_image == NULL) {
		printf("No reference image, please give a reference image (option -r)\n\n");
		print_usage(argv);
		_exit(1);
	}
	if (limit_num_layer && deflayer) {
		printf("-l and -k parameters are not compatible\n\n");
		print_usage(argv);
		_exit(1);
	}
	if (limit_num_layer && params->max_num_layer == 0) {
		printf("k should be at least 1\n");
		print_usage(argv);
		_exit(1);
	}
	if (definfolder == 0) {
		if (deflayer == 0) {
			printf("use either -l or -i option to specify layers to be computed\n\n");
			print_usage(argv);
			_exit(1);
		} else {
			if (params->nlayer == 0) {
				printf("number of layer can't be nul !\n\n");
				print_usage(argv);
				_exit(1);
			}
		}
	} else {
		if (deflayer == 1) {
			printf("choose only one of -l or -i option to specify layers to be computed\n\n");
			print_usage(argv);
			_exit(1);
		}
	}
	// printf("************* TIME MEASUREMENTS *************\n");
	fprintf(stdout, "*************** PARAMETERS **********************\n");
	if (definfolder) fprintf(stdout, "input folder :\r\t\t\t%s\n",params->infolder);
	fprintf(stdout, "reference image:\r\t\t\t%s\n",params->ref_image);
	if (params->output_dir) 
		fprintf(stdout, "output dir:\r\t\t\t%s\n",params->output_dir);
	if (params->output_mean)
		fprintf(stdout, "output mean displ:\r\t\t\tyes\n");
	else
		fprintf(stdout, "output mean displ:\r\t\t\tno\n");
	if (params->output_nacount)
		fprintf(stdout, "output NA count:\r\t\t\tyes\n");
	else
		fprintf(stdout, "output NA count:\r\t\t\tno\n");
	if (params->output_veloc)
		fprintf(stdout, "output mean velocity:\r\t\t\tyes\n");
	else
		fprintf(stdout, "output mean velocity:\r\t\t\tno\n");
	fprintf(stdout, "x resolution:\r\t\t\t%f\n",params->xres);
	fprintf(stdout, "y resolution:\r\t\t\t%f\n",params->yres);
	fprintf(stdout, "x pixel resolution:\r\t\t\t%f\n",params->px1_pas);
	fprintf(stdout, "y pixel resolution:\r\t\t\t%f\n",params->px2_pas);
#ifdef WINDOW_SIZE
	params->neighborhood = WINDOW_SIZE;
	assert((WINDOW_SIZE % 2) == 1);
#endif
	fprintf(stdout, "neighborhood:\r\t\t\t%i\n",params->neighborhood);
	fprintf(stdout, "min displacement:\r\t\t\t%f\n",params->min_displ);
	fprintf(stdout, "min correlation coef:\r\t\t\t%f\n",params->min_cc);
	fprintf(stdout, "min matches fraction:\r\t\t\t%f\n",params->min_matches);
	if (definfolder) {
		fprintf(stdout, "only multi-temporal:\r\t\t\t%s\n",params->only_multitemp?"yes":"no");
	} else {
		if (defmultitemp) printf("multi-temporal:\r\t\t\tparameter (-t) is ignored using -l parameter\n");
	}
	if (def_avail_mem) {
		printf("memory use:\r\t\t\tup to %d GiB of memory\n", (int) params->available_GiB_mem);
	} else {
		printf("memory use:\r\t\t\tparameter -g has not been set. (using default value of 2 GiB)\n");
	}
	printf("*************************************************\n");
	fprintf(stdout, "\n");

	*path_opt_ind = optind;

	if (deflayer && VERBOSE > 2) {
		fprintf(stdout, "***************** LAYERS ************************\n");
		for (int i = optind; i < argc; i++) {
			if ((i-optind)%NUM_FILES == 0) fprintf (stdout, "layer %d:\n", (i - optind)/NUM_FILES + 1);
			if ((i-optind)%NUM_FILES == 0) fprintf (stdout, "\tPX1 path:\r\t\t%s\n", argv[i]);
			if ((i-optind)%NUM_FILES == 1) fprintf (stdout, "\tPX2 path:\r\t\t%s\n", argv[i]);
			if ((i-optind)%NUM_FILES == 2) fprintf (stdout, "\tCC  path:\r\t\t%s\n", argv[i]);
			if ((i-optind)%NUM_FILES == 3) fprintf (stdout, "\tduration:\r\t\t%.2f days\n", (float)atof(argv[i]));
		}
		printf("*************************************************\n");
	}
	return 0;
}

int copy_geotiff_metadata(const char* file_dest, const char* file_src)
{
	double geo_transform[6];
	const char *projection;

	GDALAllRegister();

	// read from src
	GDALDatasetH src_dataset = GDALOpen(file_src, GA_ReadOnly);
	if(src_dataset == NULL) {
		fprintf(stderr, "Error: GDAL can't open src file %s\n", file_src);
		return 1;
	}
	if( GDALGetProjectionRef(src_dataset) != NULL ) {
		projection = GDALGetProjectionRef(src_dataset);
	} else {
		fprintf(stderr, "Error: GDAL can't read projection in %s\n", file_src);
		return 1;
	}
	if( GDALGetGeoTransform(src_dataset, geo_transform) != 0 ) {
		fprintf(stderr, "Error: GDAL can't read coordinate transformation in %s\n", file_src);
		return 1;
	}

	// write to dest
	GDALDatasetH dest_dataset = GDALOpen(file_dest, GA_Update);
	if(dest_dataset == NULL) {
		fprintf(stderr, "Error: GDAL can't open dest file %s\n", file_dest);
		return 1;
	}
	if (GDALSetProjection(dest_dataset, projection) != 0) {
		fprintf(stderr, "Error: GDAL can't set projection in %s\n", file_dest);
		return 1;
	}
	if( GDALSetGeoTransform(dest_dataset, geo_transform) != 0 ) {
		fprintf(stderr, "Error: GDAL can't set coordinate transformation in %s\n", file_dest);
		return 1;
	}

	GDALClose(src_dataset);
	GDALClose(dest_dataset);
	return 0;
}


int file_exists(char *file) 
{ 
	struct stat s; 
	if (stat(file, &s) == 0) 
		return 1 ; 
	else 
		return 0 ; 
}

int mkdir_rec(char *path)
{
	char path_dupl[2*PATH_MAX];
	char *subpath;

	if (!file_exists(path)) {
		strcpy(path_dupl, path);
		subpath = dirname(path_dupl);
		if (!file_exists(subpath)) {
			int res = mkdir_rec(subpath);
			if (res == -1) {
				fprintf(stderr, "couldn't create directory %s", subpath);
				perror(" ");
				return 1;
			} else if (res == 1) {
				return 1;
			}
		}
		return mkdir(path, 0755);
	} else {
		return 0;
	}
}

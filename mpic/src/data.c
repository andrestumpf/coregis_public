#include "data.h"

void clear_stack_f (struct stack_f *stack)
{
	memset(stack->data, 0, stack->width * stack->height * stack->nlayer * sizeof(float));
}

void clear_image_f (struct image_f *image)
{
	memset(image->data, 0, image->width * image->height * sizeof(float));
}

void clear_image_i (struct image_i *image)
{
	memset(image->data, 0, image->width * image->height * sizeof(int));
}

int allocate_domain(struct domain_info **domain, const int nb_block_x, const int nb_block_y, const int ncol, const int nrow, const int col_padding, const int row_padding)
{
	*domain = (struct domain_info *) malloc (sizeof(struct domain_info));
	if (*domain == NULL) return 1;

 	struct domain_info *dom = *domain;

	dom->nb_block_x = nb_block_x;
	dom->nb_block_y = nb_block_y;
	dom->size_block_x = (uint)ceilf((float) ncol / (float) nb_block_x);
	dom->size_block_y = (uint)ceilf((float) nrow / (float) nb_block_y);

	dom->blocks = (struct block_info *) malloc (nb_block_x * nb_block_y * sizeof(struct block_info));
	if (dom->blocks == NULL) return 1;

	for (int j = 0; j < nb_block_y; j++) {
		for (int i = 0; i < nb_block_x; i++) {
			dom->blocks[j * nb_block_x + i].num_block_x = i;
			dom->blocks[j * nb_block_x + i].num_block_y = j;

			dom->blocks[j * nb_block_x + i].size_x = (i < nb_block_x-1)? dom->size_block_x : (ncol - i * dom->size_block_x);
			dom->blocks[j * nb_block_x + i].size_y = (j < nb_block_y-1)? dom->size_block_y : (nrow - j * dom->size_block_y);

			dom->blocks[j * nb_block_x + i].first_x = (i == 0)? 1 : 0;
			dom->blocks[j * nb_block_x + i].last_x = (i == (nb_block_x - 1))? 1 : 0;

			dom->blocks[j * nb_block_x + i].first_y = (j == 0)? 1 : 0;
			dom->blocks[j * nb_block_x + i].last_y = (j == (nb_block_y - 1))? 1 : 0;

			dom->blocks[j * nb_block_x + i].x_offset_start = i * dom->size_block_x;
			dom->blocks[j * nb_block_x + i].x_offset_end = i * dom->size_block_x + dom->blocks[j * nb_block_x + i].size_x;

			dom->blocks[j * nb_block_x + i].y_offset_start = j * dom->size_block_y;
			dom->blocks[j * nb_block_x + i].y_offset_end = j * dom->size_block_y + dom->blocks[j * nb_block_x + i].size_y;
		}
	}
	return 0;
}

int get_date(struct my_date* date, const char* date_str)
{
	char year[5];
	char month[3];
	char day[3];

	strncpy(year, date_str, 4);
	year[4] = '\0';
	date->year = atoi(year);

	strncpy(month, date_str + 4, 2);
	month[2] = '\0';
	date->month = atoi(month);

	strncpy(day, date_str + 6, 2);
	day[2] = '\0';
	date->day = atoi(day);

	return 0;
}

int same_date (const struct my_date* date1, const struct my_date* date2)
{
	return (date1->day == date2->day && date1->month == date2->month && date1->year == date2->year);
}

struct stack_f *allocate_stack_f (const int nb_row, const int nb_col, const int nb_layer, const int column_padding, const int row_padding)
{
	struct stack_f *stack;
	stack = (struct stack_f*) malloc (sizeof(struct stack_f));

	stack->nrow = nb_row;
	stack->ncol = nb_col;
	stack->nlayer = nb_layer;
	stack->col_padding = column_padding;
	stack->row_padding = row_padding;

	stack->width = nb_col + 2 * column_padding;
	stack->height = nb_row + 2 * row_padding;
	stack->stride_k = stack->width * stack->height;
	stack->stride_j = stack->width * stack->nlayer;
	stack->stride_i = stack->nlayer * stack->height;

	stack->data = (float*) malloc(stack->width * stack->height * stack->nlayer * sizeof(float));

	if (stack->data == NULL) {
		fprintf(stderr, "failed to allocate stack->data !\n");
		abort();
	}

	return stack;
}

struct all_data *allocate_all_data (const int nb_row, const int nb_col, const int nb_layer, const int column_padding, const int row_padding)
{
	int i;
	struct all_data *data;
	data = (struct all_data *) malloc (sizeof(struct all_data));

	data->ew = allocate_stack_f(nb_row, nb_col, nb_layer, column_padding, row_padding);
	data->ns = allocate_stack_f(nb_row, nb_col, nb_layer, column_padding, row_padding);
	data->cc = allocate_stack_f(nb_row, nb_col, nb_layer, column_padding, row_padding);

	// no padding in images
	data->na_count = allocate_image_i(nb_row, nb_col, column_padding, row_padding);
	data->na_count_blck = allocate_image_i(nb_row, nb_col, 0, 0);
	data->mean_magn = allocate_image_f(nb_row, nb_col, 0, 0);
	data->mean_veloc_ew = allocate_image_f(nb_row, nb_col, 0, 0);
	data->mean_veloc_ns = allocate_image_f(nb_row, nb_col, 0, 0);
	data->out_mat = allocate_image_f(nb_row, nb_col, 0, 0);

	data->dates_1 = (struct my_date*) malloc(nb_layer * sizeof(struct my_date));
	data->dates_2 = (struct my_date*) malloc(nb_layer * sizeof(struct my_date));

	data->duration = (float*) malloc(nb_layer * sizeof(float));

	data->px1_paths = (char**) malloc(nb_layer * sizeof(char*));
	data->px2_paths = (char**) malloc(nb_layer * sizeof(char*));
	data->cc_paths = (char**) malloc(nb_layer * sizeof(char*));

	data->px1_is_tiled = (int*) malloc(nb_layer * sizeof(int));
	data->px2_is_tiled = (int*) malloc(nb_layer * sizeof(int));
	data->cc_is_tiled = (int*) malloc(nb_layer * sizeof(int));

	for (i = 0; i < nb_layer; i++) {
		data->px1_paths[i] = (char*) malloc(PATH_MAX * sizeof(char));
		data->px2_paths[i] = (char*) malloc(PATH_MAX * sizeof(char));
		data->cc_paths[i] = (char*) malloc(PATH_MAX * sizeof(char));
	}
	return data;
}

struct image_f *allocate_image_f(const int nb_row, const int nb_col, const int column_padding, const int row_padding)
{
	struct image_f *image;
	image = (struct image_f*) malloc (sizeof(struct image_f));

	image->nrow = nb_row;
	image->ncol = nb_col;
	image->col_padding = column_padding;
	image->row_padding = row_padding;

	image->width = nb_col + 2 * column_padding;
	image->height = nb_row + 2 * row_padding;

	image->data = (float*) malloc(image->width * image->height * sizeof(float));
	if (image->data == NULL) {
		fprintf(stderr, "failed to allocate image->data !\n");
		abort();
	}

	return image;
}

struct image_i *allocate_image_i(const int nb_row, const int nb_col, const int column_padding, const int row_padding)
{
	struct image_i *image;
	image = (struct image_i*) malloc (sizeof(struct image_i));

	image->nrow = nb_row;
	image->ncol = nb_col;
	image->col_padding = column_padding;
	image->row_padding = row_padding;

	image->width = nb_col + 2 * column_padding;
	image->height = nb_row + 2 * row_padding;

	image->data = (int*) malloc(image->width * image->height * sizeof(int));
	if (image->data == NULL) {
		fprintf(stderr, "failed to allocate image->data !\n");
		abort();
	}

	return image;
}


void free_all_data (struct all_data *data)
{
	int i, nlayer;

	nlayer = data->ew->nlayer;

	free_stack_f(data->ew);
	free_stack_f(data->ns);
	free_stack_f(data->cc);

	free_image_i(data->na_count);
	free_image_i(data->na_count_blck);
	free_image_f(data->mean_magn);
	free_image_f(data->mean_veloc_ew);
	free_image_f(data->mean_veloc_ns);
	free_image_f(data->out_mat);

	free(data->dates_1);
	free(data->dates_2);

	for (i = 0; i < nlayer; i++) {
		free(data->px1_paths[i]);
		free(data->px2_paths[i]);
		free(data->cc_paths[i]);
	}
	free(data->px1_paths);
	free(data->px2_paths);
	free(data->cc_paths);

	free(data->px1_is_tiled);
	free(data->px2_is_tiled);
	free(data->cc_is_tiled);
	free(data->duration);

	free(data);
}

void free_stack_f (struct stack_f *stack)
{
	free(stack->data);
	free(stack);
}

void free_image_f(struct image_f *image)
{
	free(image->data);
	free(image);
}

void free_image_i(struct image_i *image)
{
	free(image->data);
	free(image);
}

void free_domain(struct domain_info *domain)
{
	free(domain->blocks);
	free(domain);
}


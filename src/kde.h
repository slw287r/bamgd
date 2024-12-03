#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <cairo/cairo.h>

// Bandwidth selector
#define BANDWIDTH_X 1.0
#define BANDWIDTH_Y 1.0

typedef struct
{
	double xmin, xmax, ymin, ymax;
} dm_t;

// Bivariate KDE function
double gaussian_kernel(double u);
void biv_kde(double *x, double *y, int n, double *grid_x, double *grid_y,
		int grid_size, double *kde);
// Univariate KDE function
void univ_kde(double *data, int n, double *grid, int grid_size, double *kde,
		double bandwidth, bool horiz);
// Normalize KDE values for visualization
void norm_kde(double *kde, int size, double *min_val, double *max_val);
// Draw filled univariate KDE
void draw_filled_kde(cairo_t *cr, double *kde, int col, double *grid, int grid_size,
		int width, int height);
// Draw contours with transparency
void draw_contours(cairo_t *cr, double *kde, int col, int grid_size, int width,
		int height, int n_level);
// Round step size to a "nice" number
double calculate_tick_step(double range);
// Generalized tick drawing function
void draw_ticks(cairo_t *cr, double *grid_x, double *grid_y, int grid_size,
		int width, int height);
void draw_top_ticks(cairo_t *cr, double *grid_x, double *grid_y, int grid_size,
		int width, int height, int margin);
void draw_right_ticks(cairo_t *cr, double *grid_x, double *grid_y, int grid_size,
		int width, int height, int margin);
// Draw axes for the main bivariate plot
void draw_axes(cairo_t *cr, int width, int height);
void draw_ylab(cairo_t *cr, const char *lab, bool log, double x, double canvas_height);
void draw_legend(cairo_t *cr, int width, int margin);
void kde_plot(double *x, double *y, double *z, int n, const dm_t *dm, bool log,
		const char *png);

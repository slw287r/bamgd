#include "kde.h"
// Gaussian kernel function
double gaussian_kernel(double u)
{
	return exp(-u * u / 2.0) / sqrt(2.0 * M_PI);
}

void biv_kde(double *x, double *y, int n, double *grid_x, double *grid_y,
		int grid_size, double *kde_values)
{
	int i, j, k;
	double h_x = BANDWIDTH_X, h_y = BANDWIDTH_Y;
	for (i = 0; i < grid_size; ++i)
	{
		for (j = 0; j < grid_size; ++j)
		{
			double gx = grid_x[i], gy = grid_y[j], sum = 0.0f;
			for (k = 0; k < n; ++k)
			{
				double u = (gx - x[k]) / h_x, v = (gy - y[k]) / h_y;
				sum += gaussian_kernel(u) * gaussian_kernel(v);
			}
			//kde_values[i * grid_size + (grid_size - j)] = sum / (n * h_x * h_y);
			kde_values[i * grid_size + j] = sum / (n * h_x * h_y);
		}
	}
	/* debug kde
	//for (i = 0; i < grid_size; ++i)
	for (i = 0; i < 1; ++i)
		for (j = 0; j < grid_size; ++j)
			printf("%d\t%d\t%f\n", i, j, kde_values[i * grid_size + (grid_size - j)]);
	*/
}

void univ_kde(double *data, int n, double *grid, int grid_size, double *kde_values,
		double bandwidth, bool horiz)
{
	int i, k;
	for (i = 0; i < grid_size; ++i)
	{
		double gx = grid[i];
		double sum = 0.0;
		for (k = 0; k < n; ++k)
		{
			double u = (gx - data[k]) / bandwidth;
			sum += gaussian_kernel(u);
		}
		kde_values[horiz ? i : grid_size - i] = sum / (n * bandwidth);
	}
}

void norm_kde_values(double *kde_values, int size, double *min_val, double *max_val)
{
	int i;
	*min_val = *kde_values;
	*max_val = *kde_values;
	for (i = 0; i < size; ++i)
	{
		*min_val = fmin(*min_val, kde_values[i]);
		*max_val = fmax(*max_val, kde_values[i]);
	}
	for (i = 0; i < size; ++i)
		kde_values[i] = (kde_values[i] - *min_val) / (*max_val - *min_val);
}

void draw_filled_kde(cairo_t *cr, double *kde_values, double *grid, int grid_size,
		int width, int height)
{
	int i;
	cairo_set_source_rgba(cr, 199 / 255.0, 219 / 255.0, 217 / 255.0, 0.6); // Steelblue fill with transparency
	cairo_move_to(cr, 0, height); // Start at the base of the plot
	for (i = 0; i < grid_size; ++i)
	{
		double gx = (double)i / (grid_size - 1) * width; // Scaled grid position
		double gy = kde_values[i] * height;       // Scaled density value
		cairo_line_to(cr, gx, height - gy * 0.9);
	}
	cairo_line_to(cr, width, height);
	cairo_close_path(cr);
	cairo_fill(cr);
	cairo_set_source_rgba(cr, 101 / 255.0, 175 / 255.0, 144 / 255.0, 0.9); // Steelblue line with transparency
	cairo_move_to(cr, 0, (1 - kde_values[1] * 0.9) * height);
	for (i = 1; i < grid_size; ++i)
	{
		double gx = (double)i / (grid_size - 1) * width;
		double gy = kde_values[i] * height;
		cairo_line_to(cr, gx, height - gy * 0.9);
	}
	cairo_stroke(cr);
}

// Draw contours with transparency
void draw_contours(cairo_t *cr, double *kde_values, int grid_size, int width,
		int height, int n_level)
{
	int i, j, level;
	double cell_width = (double)width / grid_size;
	double cell_height = (double)height / grid_size;
	for (level = 1; level <= n_level; ++level)
	{
		double threshold = (double)level / n_level;
		cairo_set_source_rgba(cr, 9 / 255.0, 15 / 255.0, 13 / 255.0, 0.1 * level); // Steelblue color with varying alpha
		for (i = 0; i < grid_size; ++i)
		{
			for (j = 0; j < grid_size; ++j)
			{
				double value = kde_values[i * grid_size + j];
				if (value >= threshold)
				{
					cairo_rectangle(cr, i * cell_width, (grid_size - j) * cell_height, cell_width,
							cell_height);
					cairo_fill(cr);
				}
			}
		}
	}
}

// Round step size to a "nice" number
double calculate_tick_step(double range)
{
	double base_step = pow(10, floor(log10(range))); // Base step (e.g., 1, 10, 0.1)
	double relative_step = range / base_step;
	// Adjust step size to "nice" intervals
	if (relative_step < 2)
		return base_step / 5.0;  // Fine steps
	if (relative_step < 5)
		return base_step / 2.0;  // Medium steps
	return base_step;
}

// Generalized tick drawing function
void draw_ticks(cairo_t *cr, double *grid_x, double *grid_y, int grid_size,
		int width, int height)
{
	double x, y;
	cairo_text_extents_t ext;
	cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0); // Black for ticks and labels
	cairo_set_line_width(cr, 1.5);
	cairo_set_font_size(cr, 18.0);
	// Compute ranges and dynamic tick steps
	double x_min = grid_x[0], x_max = grid_x[grid_size - 1];
	double y_min = grid_y[0], y_max = grid_y[grid_size - 1];
	double x_range = x_max - x_min;
	double y_range = y_max - y_min;
	double x_step = calculate_tick_step(x_range);
	double y_step = calculate_tick_step(y_range);
	// X-axis ticks
	for (x = ceil(x_min / x_step) * x_step; x <= x_max; x += x_step)
	{
		double x_pos = (x - x_min) / x_range * width;
		cairo_move_to(cr, x_pos, height);
		cairo_set_source_rgb(cr, 1, 1, 1);
		cairo_line_to(cr, x_pos, 0);
		cairo_stroke(cr);
		// Add tick label
		char label[NAME_MAX];
		snprintf(label, sizeof(label), "%g", x);
		cairo_text_extents(cr, label, &ext);
		cairo_set_source_rgb(cr, 0, 0, 0);
		cairo_move_to(cr, x_pos - ext.width / 2 - ext.x_bearing, height + ext.height * 2); // Adjust position for better alignment
		cairo_show_text(cr, label);
	}
	// Y-axis ticks
	// get precision
	int prec = 0;
	for (y = ceil(y_min / y_step) * y_step; y <= y_max; y += y_step)
	{
		char label[NAME_MAX], *p;
		snprintf(label, sizeof(label), "%g", y);
		prec = fmax(prec, (p = strchr(label, '.')) ? strlen(p + 1) : 0);
	}
	for (y = ceil(y_min / y_step) * y_step; y <= y_max; y += y_step)
	{
		double y_pos = (y - y_min) / y_range * height;
		cairo_move_to(cr, 0, height - y_pos);
		cairo_set_source_rgb(cr, 1, 1, 1);
		cairo_line_to(cr, width, height - y_pos);
		cairo_stroke(cr);
		// Add tick label
		char label[NAME_MAX];
		snprintf(label, sizeof(label), "%.*f", prec, y);
		cairo_set_source_rgb(cr, 0, 0, 0);
		cairo_text_extents(cr, label, &ext);
		cairo_move_to(cr, -ext.width -ext.x_bearing - ext.height / 2, height - y_pos + ext.height / 2); // Adjust position for better alignment
		cairo_show_text(cr, label);
	}
}

void draw_top_ticks(cairo_t *cr, double *grid_x, double *grid_y, int grid_size,
		int width, int height, int margin)
{
	double x;
	cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0); // Black for ticks and labels
	cairo_set_line_width(cr, 1.5);
	// Compute ranges and dynamic tick steps
	double x_min = grid_x[0], x_max = grid_x[grid_size - 1];
	double x_range = x_max - x_min;
	double x_step = calculate_tick_step(x_range);
	// X-axis ticks
	for (x = ceil(x_min / x_step) * x_step; x <= x_max; x += x_step)
	{
		double x_pos = (x - x_min) / x_range * width;
		cairo_move_to(cr, x_pos, margin);
		cairo_set_source_rgb(cr, 1, 1, 1);
		cairo_line_to(cr, x_pos, 0);
		cairo_stroke(cr);
	}
}

void draw_right_ticks(cairo_t *cr, double *grid_x, double *grid_y, int grid_size,
		int width, int height, int margin)
{
	double y;
	cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0); // Black for ticks and labels
	cairo_set_line_width(cr, 1.5);
	// Compute ranges and dynamic tick steps
	double y_min = grid_y[0], y_max = grid_y[grid_size - 1];
	double y_range = y_max - y_min;
	double y_step = calculate_tick_step(y_range);
	// Y-axis ticks
	for (y = ceil(y_min / y_step) * y_step; y <= y_max; y += y_step)
	{
		double y_pos = (y - y_min) / y_range * height;
		cairo_move_to(cr, width + margin, height - y_pos);
		cairo_set_source_rgb(cr, 1, 1, 1);
		cairo_line_to(cr, width + margin * 2, height - y_pos);
		cairo_stroke(cr);
	}
}

// Draw axes for the main bivariate plot
void draw_axes(cairo_t *cr, int width, int height)
{
	cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0); // Black axis lines
	cairo_set_line_width(cr, 1.5);
	// X-axis line
	cairo_move_to(cr, 0, height); // Start at bottom-left corner
	cairo_line_to(cr, width, height); // Draw to bottom-right corner
	cairo_stroke(cr);
	// Y-axis line
	cairo_move_to(cr, 0, 0); // Start at top-left corner
	cairo_line_to(cr, 0, height); // Draw to bottom-left corner
	cairo_stroke(cr);
}

void draw_ylab(cairo_t *cr, const char *lab, double x, double canvas_height)
{
	cairo_set_font_size(cr, 20.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_text_extents_t ext;
	cairo_text_extents(cr, lab, &ext);
	double y_center = canvas_height / 2.0;
	double x_pos = x - (ext.height / 2.0);
	double y_pos = y_center + (ext.width / 2.0);
	cairo_save(cr);
	cairo_translate(cr, x_pos, y_pos);
	cairo_rotate(cr, -M_PI / 2.0);
	cairo_move_to(cr, 0, 0);
	cairo_show_text(cr, lab);
	cairo_restore(cr);
}

// Main program
void kde_plot(double *x, double *y, int n, const dm_t *dm, bool log, const char *png)
{
	int i;
	// Grid definition
	int points = 512;
	int grid_size = points;
	double grid_x[grid_size];
	double grid_y[grid_size];
	double kde_values[grid_size * grid_size];
	double kde_x[grid_size], kde_y[grid_size];
	// Create grid points
	for (i = 0; i < grid_size; ++i)
	{
		grid_x[i] = dm->xmin + i * (dm->xmax - dm->xmin) / (grid_size - 1);
		grid_y[i] = dm->ymin + i * (dm->ymax - dm->ymin) / (grid_size - 1);
	}
	// Compute KDE
	biv_kde(x, y, n, grid_x, grid_y, grid_size, kde_values);
	univ_kde(x, n, grid_x, grid_size, kde_x, BANDWIDTH_X, true);
	univ_kde(y, n, grid_y, grid_size, kde_y, BANDWIDTH_Y, false);
	// Normalize for visualization
	double min_val, max_val;
	norm_kde_values(kde_values, grid_size * grid_size, &min_val, &max_val);
	norm_kde_values(kde_x, grid_size, &min_val, &max_val);
	norm_kde_values(kde_y, grid_size, &min_val, &max_val);
	// Visualization canvas
	double width = points, height = points;
	double margin = 64, n_level = 10;
	cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32,
			width + 2 * margin, height + 2 * margin);
	cairo_t *cr = cairo_create(surface);
	// white background
	cairo_save(cr);
	cairo_set_source_rgb(cr, 1, 1, 1);
	cairo_paint(cr);
	cairo_restore(cr);
	// Draw bivariate KDE with contours
	cairo_save(cr);
	cairo_translate(cr, margin, margin); // Offset for main plot
	cairo_rectangle(cr, 0, 0, width, height); // background
	cairo_set_source_rgb(cr, 234 / 255.0, 234 / 255.0, 241 / 255.0);
	cairo_fill(cr);
	draw_ticks(cr, grid_x, grid_y, grid_size, width, height);
	// draw_axes(cr, width, height);
	draw_contours(cr, kde_values, grid_size, width, height, n_level);
	cairo_restore(cr);
	// Draw top univariate KDE plot
	cairo_save(cr);
	cairo_translate(cr, margin, 0); // Position top plot above main plot
	cairo_rectangle(cr, 0, margin * 0.05, width, margin * 0.9); // background
	cairo_set_source_rgb(cr, 234 / 255.0, 234 / 255.0, 241 / 255.0);
	cairo_fill(cr);
	draw_top_ticks(cr, grid_x, grid_y, grid_size, width, height, margin);
	draw_filled_kde(cr, kde_x, grid_x, grid_size, width, margin * 0.95);
	cairo_restore(cr);
	// Draw right univariate KDE plot
	cairo_save(cr);
	cairo_translate(cr, width + margin * 2, margin); // Position right plot beside main plot
	cairo_rotate(cr, M_PI / 2);
	cairo_rectangle(cr, 0, margin * 0.05, height, margin * 0.9); // background
	cairo_set_source_rgb(cr, 234 / 255.0, 234 / 255.0, 241 / 255.0);
	cairo_fill(cr);
	cairo_restore(cr);
	// right ticks
	cairo_save(cr);
	cairo_translate(cr, 0, margin); // Offset for main plot
	draw_right_ticks(cr, grid_x, grid_y, grid_size, width, height, margin * 0.95);
	cairo_restore(cr);
	// density
	cairo_save(cr);
	cairo_translate(cr, width + margin * 2, margin); // Position right plot beside main plot
	cairo_rotate(cr, M_PI / 2);
	draw_filled_kde(cr, kde_y, grid_y, grid_size, height, margin * 0.95);
	cairo_restore(cr);
	// axis labels
	cairo_text_extents_t ext;
	char *xlab = "GC (%)", *ylab = log ? "Depth (log₁₀×)" : "Depth (×)";
	cairo_set_font_size(cr, 20.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_text_extents(cr, xlab, &ext);
	cairo_move_to(cr, width / 2.0 + margin - (ext.width / 2.0 + ext.x_bearing),
			height + margin * 2 - ext.height);
	cairo_show_text(cr, xlab);
	cairo_stroke(cr);
	draw_ylab(cr, ylab, margin / 2.0, height + margin * 2);
	cairo_surface_write_to_png(surface, png);
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

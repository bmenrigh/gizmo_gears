#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <quadmath.h>
#include <assert.h>

#include <png.h>

#define ORD_LIMIT 50000000


struct samples {
    uint32_t count;
    double order;
};

struct render_ctx {
    int n;
    double r;
    int img_w, img_h;
    double xmin, xmax, ymin, ymax;
    double pwidth, pheight, half_pwidth, half_pheight, pradius;
    __complex128 rot_a, rot_b;
    __float128 epsilon;
    struct samples *grid;
};


void write_png_file(char *filename, int width, int height, png_bytep image_data) {

  FILE *fp = fopen(filename, "wb");
  assert(fp != NULL);

  png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  assert(png != NULL);

  png_infop info = png_create_info_struct(png);
  assert(info != NULL);

  /* PNG "exception handling" */
  assert(setjmp(png_jmpbuf(png)) == 0);

  png_init_io(png, fp);

  png_set_IHDR(
    png,
    info,
    width, height,
    8,
    PNG_COLOR_TYPE_RGB,
    PNG_INTERLACE_NONE,
    PNG_COMPRESSION_TYPE_DEFAULT,
    PNG_FILTER_TYPE_DEFAULT
  );
  png_write_info(png, info);

  assert(image_data != NULL);

  for (int y = 0; y < height; y++) {
      png_write_row(png, &(image_data[y * (width * 3)]));
  }

  png_write_end(png, NULL);

  fclose(fp);

  png_destroy_write_struct(&png, &info);
}


void ctx_to_png(struct render_ctx *ctx, char *name) {

    png_bytep image_data = calloc(ctx->img_h * ctx->img_w * 3, sizeof(png_byte));

    /* Find min and max order */
    double log_max_order = 0;
    double log_min_order = ORD_LIMIT;
    for (int y = 0; y < ctx->img_h; y++) {
        for (int x = 0; x < ctx->img_w; x++) {
            int o = y * ctx->img_w + x;
            if (ctx->grid[o].count > 0) {
                double log_avg_order = fabs((double)ctx->grid[o].order / (double)ctx->grid[o].count);
                if (log_avg_order > log_max_order) {
                    log_max_order = log_avg_order;
                }
                if (log_avg_order < log_min_order) {
                    log_min_order = log_avg_order;
                }
            }
        }
    }

    double max_order = exp(log_max_order);
    double min_order = exp(log_min_order);

    fprintf(stderr, "log min: %f; min: %f\n", log_min_order, min_order);
    fprintf(stderr, "log max: %f; max: %f\n", log_max_order, max_order);

    /* color pixles scaled by max order */
    for (int y = 0; y < ctx->img_h; y++) {
        for (int x = 0; x < ctx->img_w; x++) {
            int o = y * ctx->img_w + x;
            if (ctx->grid[o].count > 0) {
                int neg = 0;
                double log_avg_order = (double)ctx->grid[o].order / (double)ctx->grid[o].count;
                if (log_avg_order < 0) {
                    neg = 1;
                    log_avg_order = fabs(log_avg_order);
                }

                double avg_order = exp(log_avg_order);

                double v = atan2(avg_order - min_order, 1.0) / atan2(max_order - min_order, 1.0);

                /* int grey = (int)floor(((avg_order - min_order) / (max_order - min_order)) * 255.0);*/
                /*int grey = (int)floor(v * 255);


                if (neg == 0) {
                    image_data[o * 3 + 0] = grey;
                } else {
                    image_data[o * 3 + 1] = grey;
                }
                image_data[o * 3 + 2] = grey;*/


                /* delta colors */
                v = v * M_PI_2;

                if (neg == 0) {
                    image_data[o * 3 + 0] = (int)round(sin(v) * 255.0);
                    image_data[o * 3 + 1] = (int)round((1.0 - sin(v * 2.0)) * 255.0);
                    image_data[o * 3 + 2] = (int)round(cos(v) * 255.0);
                } else {
                    image_data[o * 3 + 0] = (int)round((1.0 - sin(v)) * 255.0);
                    image_data[o * 3 + 1] = (int)round(sin(v * 2.0) * 255.0);
                    image_data[o * 3 + 2] = (int)round((1.0 - cos(v)) * 255.0);
                }
            }
        }
    }

    write_png_file(name, ctx->img_w, ctx->img_h, image_data);

    free(image_data);
}


int point_to_xy(struct render_ctx *ctx, __complex128 p, int *x, int *y) {

    double px = (double)(__real__ p);
    double py = (double)(__imag__ p);

    /*fprintf(stderr, "point (x, y): (%.5f, %.5f)\n", px, py);*/

    if ((px + ctx->half_pwidth < ctx->xmin) ||
        (px - ctx->half_pwidth > ctx->xmax) ||
        (py + ctx->half_pheight < ctx->ymin) ||
        (py - ctx->half_pheight > ctx->ymax)) {
        return -1;
    }

    *x = (int)floor(((px - ctx->xmin) / ctx->pwidth) + 0.5);
    *y = (int)floor(((py - ctx->ymin) / ctx->pheight) + 0.5);

    return 0;
}


__complex128 point_from_xy(struct render_ctx *ctx, int x, int y) {

    __complex128 p;

    __real__ p = (ctx->xmin + (ctx->pwidth * x));
    __imag__ p = (ctx->ymin + (ctx->pheight * y));

    return p;
}


__complex128 point_rand_offset(struct render_ctx *ctx) {

    __complex128 p;
    __real__ p = (__float128)((((double)rand() / (double)0x80000000) - 0.5) * ctx->pwidth);
    __imag__ p = (__float128)((((double)rand() / (double)0x80000000) - 0.5) * ctx->pheight);

    return p;
}


__complex128 point_from_xy_rand(struct render_ctx *ctx, int x, int y) {

    return point_from_xy(ctx, x, y) + point_rand_offset(ctx);
}


int point_equal_double(__complex128 p, __complex128 q) {

    if (((double)(__real__ p) == (double)(__real__ q)) &&
        ((double)(__imag__ p) == (double)(__imag__ q))) {

        return 1;
    } else {
        return 0;
    }
}


int point_equal_epsilon(struct render_ctx *ctx, __complex128 p, __complex128 q) {

    /* if (cabsq(p - q) < ctx->epsilon) { */
    if (fabsq(__real__ p - __real__ q) + fabsq(__imag__ p - __imag__ q)  < ctx->epsilon) {
        return 1;
    } else {
        return 0;
    }
}


int point_in_a(struct render_ctx *ctx, __complex128 p) {

    __complex128 np = p;

    __real__ np += (__float128)1;

    if ((double)cabsq(np) <= ctx->r) {
        return 1;
    } else {
        return 0;
    }
}


int point_in_b(struct render_ctx *ctx, __complex128 p) {

    __complex128 np = p;

    __real__ np -= (__float128)1;

    if ((double)cabsq(np) <= ctx->r) {
        return 1;
    } else {
        return 0;
    }
}


int xy_on_border(struct render_ctx *ctx, int x, int y) {

    __complex128 p = point_from_xy(ctx, x, y);
    __complex128 np;
    double dist;

    /* Check distace to A */
    np = p;
    __real__ np += (__float128)1;

    dist = (double)cabsq(np);
    if ((dist <= ctx->r + ctx->pradius) && (dist >= ctx->r - ctx->pradius)) {
        return 1;
    }

    /* Check distace to B */
    np = p;
    __real__ np -= (__float128)1;

    dist = (double)cabsq(np);
    if ((dist <= ctx->r + ctx->pradius) && (dist >= ctx->r - ctx->pradius)) {
        return 1;
    }

    return 0;
}


int point_in_puzzle(struct render_ctx *ctx, __complex128 p) {

    if ((point_in_a(ctx, p) == 1) || (point_in_b(ctx, p) == 1)) {
        return 1;
    } else {
        return 0;
    }
}


__complex128 turn_a(struct render_ctx *ctx, __complex128 p) {

    __complex128 np = p;
    __real__ np += (__float128)1;

    np *= ctx->rot_a;

    __real__ np -= (__float128)1;

    return np;
}


__complex128 turn_b(struct render_ctx *ctx, __complex128 p) {

    __complex128 np = p;
    __real__ np -= (__float128)1;

    np *= ctx->rot_b;

    __real__ np += (__float128)1;

    return np;
}


double point_order(struct render_ctx *ctx, __complex128 p, int32_t limit, uint8_t *visited, uint8_t *visited_m) {

    if (point_in_puzzle(ctx, p) != 1) {
        return 0;
    }

    int x, y;
    __complex128 op = p; /* Original p */
    uint8_t step = 0;
    int32_t count = 0;
    int32_t count_a = 0;
    int32_t count_b = 0;
    do {
        if (step == 0) {
            if (point_in_a(ctx, p) == 1) {
                p = turn_a(ctx, p);
                count_a++;
            }
        } else {
            if (point_in_b(ctx, p) == 1) {
                p = turn_b(ctx, p);
                count_b++;
            }
        }

        if (step == 0) {
            /* Only track after a0123456789101112 is done */
            assert(point_to_xy(ctx, p, &x, &y) == 0);
            visited[y * ctx->img_w + x] = 1;

            __complex128 p_m;
            __real__ p_m = -1.0Q * __real__ p;
            __imag__ p_m = -1.0Q * __imag__ p;

            assert(point_to_xy(ctx, p_m, &x, &y) == 0);
            visited_m[y * ctx->img_w + x] = 1;
        }

        step ^= 1; /* toggle between a and b */
        count++;

        if (count > limit) {
            return NAN;
        }

        if (count % 10000000 == 0) {
            fprintf(stderr, "Done %d turns\n", count);
        }
    } while ((count < ctx->n) || (point_equal_epsilon(ctx, op, p) != 1));

    /* regular order */
    /*return count;*/

    /* delta order a/b stuff */
    if ((count_a > 0) && (count_b > 0)) {
        return (log((double)count_a) - log((double)count_b)); /* log(a/b) */
    } else {
        return NAN;
    }
}


void point_sample(struct render_ctx *ctx, __complex128 p, int32_t limit, uint8_t *visited, uint8_t *visited_m) {

    if (point_in_puzzle(ctx, p) != 1) {
        return;
    }

    memset(visited, 0, ctx->img_w * ctx->img_h * sizeof(uint8_t));
    memset(visited_m, 0, ctx->img_w * ctx->img_h * sizeof(uint8_t));
    double ord = point_order(ctx, p, limit, visited, visited_m);

    if (isnan(ord) == 0) {
        for (int y = 0; y < ctx->img_h; y++) {
            for (int x = 0; x < ctx->img_w; x++) {
                int o = y * ctx->img_w + x;
                if (visited[o] == 1) {
                    ctx->grid[o].count += 1;
                    ctx->grid[o].order += ord;
                }

                if (visited_m[o] == 1) {
                    ctx->grid[o].count += 1;
                    ctx->grid[o].order -= ord;
                }
            }
        }
    } else {
        /* Fill in failed pixels sampled with limit */

        /* for (int y = 0; y < ctx->img_h; y++) { */
        /*     for (int x = 0; x < ctx->img_w; x++) { */
        /*         if (visited[y * ctx->img_w + x] == 1) { */
        /*             ctx->grid[y * ctx->img_w + x].count += 1; */
        /*             ctx->grid[y * ctx->img_w + x].order += log(limit); */
        /*         } */
        /*     } */
        /* } */

    }
}


void xy_sample(struct render_ctx *ctx, int x, int y, uint32_t n, uint32_t m, int32_t limit, uint8_t *visited, uint8_t *visited_m) {


    int o = y * ctx->img_w + x;
    for (uint32_t i = 0; i < n; i++) {
        if (ctx->grid[o].count < m) {
            __complex128 p = point_from_xy_rand(ctx, x, y);

            point_sample(ctx, p, limit, visited, visited_m);
        } else {
            break;
        }
    }
}


int main (void) {


    struct render_ctx *ctx = calloc(1, sizeof(struct render_ctx));
    ctx->n = 12;
    ctx->r = sqrt(2.0);
    ctx->epsilon = 1e-16Q;

    double goalw = 10240;
    double scalef = goalw / ((2.0 * ctx->r) + 2.0);
    ctx->img_w = (int)floor(((2.0 * ctx->r) + 2.0) * scalef);
    ctx->img_h = (int)floor(2.0 * ctx->r * scalef);
    ctx->xmin = -1.0 - ctx->r;
    ctx->xmax = 1.0 + ctx->r;
    ctx->ymin = -1.0 * ctx->r;
    ctx->ymax = 1.0 * ctx->r;
    ctx->pwidth = (ctx->xmax - ctx->xmin) / ((double)ctx->img_w - 1.0);
    ctx->pheight = (ctx->ymax - ctx->ymin) / ((double)ctx->img_h - 1.0);
    ctx->half_pwidth = ctx->pwidth / 2.0;
    ctx->half_pheight = ctx->pheight / 2.0;
    ctx->grid = calloc(ctx->img_w * ctx->img_h, sizeof(struct samples));
    ctx->pradius = sqrt(pow(ctx->half_pwidth, 2.0) + pow(ctx->half_pheight, 2.0));

    __real__ ctx->rot_a = cosq(M_PIq * ((__float128)2 / (__float128)ctx->n));
    __imag__ ctx->rot_a = sinq(M_PIq * ((__float128)2 / (__float128)ctx->n));
    ctx->rot_b = conjq(ctx->rot_a);

    fprintf(stderr, "Image size %dx%d\n", ctx->img_w, ctx->img_h);

    fprintf(stderr, "Image box (%.5f, %.5f) to (%.5f, %.5f)\n",
            ctx->xmin, ctx->ymin,
            ctx->xmax, ctx->ymax);


    /* Test point -> xy and xy -> point */
    int tx = 13, ty = 31;
    int nx, ny;
    __complex128 tp;
    for (int i = 0; i < 1000; i++) {
        nx = 0; ny = 0;
        tp = point_from_xy_rand(ctx, tx, ty);
        assert(point_to_xy(ctx, tp, &nx, &ny) == 0);

        if ((nx != tx) || (ny != ty)) {
            fprintf(stderr, "Failed for point (%.5f, %.5f)\n", (double)__real__ tp, (double)__imag__ tp);
        }
    }


    /*__complex128 p;
    __real__ p = -0.000001;
    __imag__ p = -0.000001;*/
    /*int x, y;
    assert(point_to_xy(ctx, p, &x, &y) == 0);
    fprintf(stderr, "(x, y): (%d, %d)\n", x, y);
    __complex128 q = point_from_xy(ctx, x, y);
    fprintf(stderr, "(x, y): (%.5f, %.5f)\n", (double)__real__ q, (double)__imag__ q);*/

    uint8_t *visited = malloc(ctx->img_w * ctx->img_h * sizeof(uint8_t));
    uint8_t *visited_m = malloc(ctx->img_w * ctx->img_h * sizeof(uint8_t));

    fprintf(stderr, "== SAMPLING BORDER PIXELS ==\n");
    for (int y = 0; y < ctx->img_h; y++) {
        if (y % 10 == 0) {
            fprintf(stderr, "Working on row %d of %d\n", y, ctx->img_h);
        }
        for (int x = 0; x < ctx->img_w; x++) {

            if (xy_on_border(ctx, x, y) == 1) {
                xy_sample(ctx, x, y, 64, 4, ORD_LIMIT, visited, visited_m);
            }
        }
    }

    fprintf(stderr, "== SAMPLING DISC PIXELS ==\n");
    for (int y = 0; y < ctx->img_h; y++) {
        if (y % 10 == 0) {
            fprintf(stderr, "Working on row %d of %d\n", y, ctx->img_h);
        }
        for (int x = 0; x < ctx->img_w; x++) {

            xy_sample(ctx, x, y, 64, 8, ORD_LIMIT, visited, visited_m);
        }
    }


    free(visited);
    free(visited_m);

    ctx_to_png(ctx, "/tmp/test.png");


    free(ctx->grid);
    free(ctx);

    return 0;
}

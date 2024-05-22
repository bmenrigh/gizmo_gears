#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <quadmath.h>
#include <assert.h>
#include <setjmp.h>
#include <pthread.h>

#include <png.h>


#define ORD_LIMIT (50 * 1000 * 1000)
#define NUM_THREADS 24
#define LOG_SCALE 0x10000
#define GROW_VISITED 1024


struct thread_ctx {
    pthread_t tid;
    int tnum;
    int num_threads;
    struct render_ctx *ctx;
};

struct samples {
    uint32_t count;
    int64_t scaled_log_order;
    /*uint64_t ord_a, ord_b;*/
};

struct visited_ctx {
    uint32_t limit;
    uint8_t *visited, *visited_m;
    uint32_t *vx, *vy, *vx_m, *vy_m;
    int32_t vsize, vused, vsize_m, vused_m;
    uint32_t ord_a, ord_b;
};

struct render_ctx {
    uint32_t n;
    double r;
    double r_sq;
    int img_w, img_h;
    double xmin, xmax, ymin, ymax;
    double pwidth, pheight, half_pwidth, half_pheight, pradius;
    __complex128 rot_a, rot_b;
    __float128 epsilon;
    struct samples *grid;
    pthread_mutex_t *grid_mutex;
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

    /* Check triangle inequality first */
    if ((double)(fabsq(__real__ np) + fabsq(__imag__ np)) < ctx->r) {
        return 1;

        /* Else check squared pythagorean */
    } else if ((double)(((__real__ np) * (__real__ np)) +
                        ((__imag__ np) * (__imag__ np))) < ctx->r_sq) {
        return 1;
    } else {
        return 0;
    }
}


int point_in_b(struct render_ctx *ctx, __complex128 p) {

    __complex128 np = p;
    __real__ np -= (__float128)1;

    /* Check triangle inequality first */
    if ((double)(fabsq(__real__ np) + fabsq(__imag__ np)) < ctx->r) {
        return 1;

        /* Else check squared pythagorean */
    } else if ((double)(((__real__ np) * (__real__ np)) +
                        ((__imag__ np) * (__imag__ np))) < ctx->r_sq) {
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


void add_visited_xy(struct render_ctx *ctx, struct visited_ctx *vctx, int x, int y, int x_m, int y_m) {

    uint32_t o = y * ctx->img_w + x;
    if (vctx->visited[o] == 0) {
        vctx->visited[o] = 1;

        /* grow x/y array */
        if (vctx->vused >= (vctx->vsize - 1)) {
            vctx->vx = realloc(vctx->vx, (vctx->vsize + GROW_VISITED) * sizeof(uint32_t));
            assert(vctx->vx != NULL);

            vctx->vy = realloc(vctx->vy, (vctx->vsize + GROW_VISITED) * sizeof(uint32_t));
            assert(vctx->vy != NULL);

            vctx->vsize += GROW_VISITED;
        }

        vctx->vx[vctx->vused] = x;
        vctx->vy[vctx->vused] = y;
        vctx->vused += 1;
    } else if (vctx->visited[o] < 32) {
        /* Let visited counter grow as high as 32 per pixel per point sampled */
        vctx->visited[o] =+ 1;
    }

    uint32_t o_m = y_m * ctx->img_w + x_m;
    if (vctx->visited_m[o_m] == 0) {
        vctx->visited_m[o_m] = 1;

        /* grow x/y array */
        if (vctx->vused_m >= (vctx->vsize_m - 1)) {
            vctx->vx_m = realloc(vctx->vx_m, (vctx->vsize_m + GROW_VISITED) * sizeof(uint32_t));
            assert(vctx->vx_m != NULL);

            vctx->vy_m = realloc(vctx->vy_m, (vctx->vsize_m + GROW_VISITED) * sizeof(uint32_t));
            assert(vctx->vy_m != NULL);

            vctx->vsize_m += GROW_VISITED;
        }

        vctx->vx_m[vctx->vused_m] = x_m;
        vctx->vy_m[vctx->vused_m] = y_m;
        vctx->vused_m += 1;
    } else if (vctx->visited_m[o_m] < 32) {
        /* Let visited counter grow as high as 32 per pixel per point sampled */
        vctx->visited_m[o_m] =+ 1;
    }
}


double point_order(struct render_ctx *ctx, __complex128 p, struct visited_ctx *vctx) {

    if (point_in_puzzle(ctx, p) != 1) {
        return 0;
    }

    int x, y, x_m, y_m;
    __complex128 op = p; /* Original p */
    uint8_t step = 0;
    uint32_t count = 0;
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

            __complex128 p_m;
            __real__ p_m = -1.0Q * __real__ p;
            __imag__ p_m = -1.0Q * __imag__ p;

            assert(point_to_xy(ctx, p_m, &x_m, &y_m) == 0);

            add_visited_xy(ctx, vctx, x, y, x_m, y_m);
        }

        step ^= 1; /* toggle between a and b */
        count++;

        if (count > vctx->limit) {

            /* Try to salvage this point if it's extremely close to 0 */
            /* Within one loop around of each other */
            if (abs(count_a - count_b) <= ctx->n) {
                fprintf(stderr, "Salvaged point at limit\n");
                return 0.0;
            } else if (count > ORD_LIMIT) {
                return 0.0; /* Just assume they were equal as a speedup hack */
            } else {
                /* We didn't try long enough to be sure */
                return NAN;
            }

        }

        if (count % 10000000 == 0) {
            fprintf(stderr, "Done %d turns\n", count);
        }
    } while ((count < ctx->n) || (point_equal_epsilon(ctx, op, p) != 1));

    /* regular order */
    /*return count;*/

    /* delta order a/b stuff */
    if ((count_a > 0) && (count_b > 0)) {

        /* *ord_a = count_a; */
        /* *ord_b = count_b; */

        return (log((double)count_a) - log((double)count_b)); /* log(a/b) */
    } else {
        return NAN;
    }
}


void point_sample(struct render_ctx *ctx, __complex128 p, struct visited_ctx *vctx) {

    if (point_in_puzzle(ctx, p) != 1) {
        return;
    }

    /* memset(vctx->visited, 0, ctx->img_w * ctx->img_h * sizeof(uint8_t)); */
    /* memset(vctx->visited_m, 0, ctx->img_w * ctx->img_h * sizeof(uint8_t)); */
    vctx->vused = 0;
    vctx->vused_m = 0;

    double ord = point_order(ctx, p, vctx);

    if (isnan(ord) == 0) {

        int64_t scaled_ord = (int64_t)round(ord * (double)LOG_SCALE);

        for (int i = 0; i < vctx->vused; i++) {
            int o = vctx->vy[i] * ctx->img_w + vctx->vx[i];

            __sync_add_and_fetch(&(ctx->grid[o].count), vctx->visited[o]);
            /* __sync_add_and_fetch(&(ctx->grid[o].ord_a), ord_a); */
            /* __sync_add_and_fetch(&(ctx->grid[o].ord_b), ord_b); */
            __sync_add_and_fetch(&(ctx->grid[o].scaled_log_order), scaled_ord * vctx->visited[o]);

            vctx->visited[o] = 0; /* clear this visit */
        }

        for (int i = 0; i < vctx->vused_m; i++) {
            int o_m = vctx->vy_m[i] * ctx->img_w + vctx->vx_m[i];

            __sync_add_and_fetch(&(ctx->grid[o_m].count), vctx->visited_m[o_m]);
            /* __sync_add_and_fetch(&(ctx->grid[o].ord_a), ord_b); /\* a/b swapped *\/ */
            /* __sync_add_and_fetch(&(ctx->grid[o].ord_b), ord_a); */
            __sync_sub_and_fetch(&(ctx->grid[o_m].scaled_log_order), scaled_ord * vctx->visited_m[o_m]);

            vctx->visited_m[o_m] = 0; /* clear this visit */
        }

    } else {
        /* Gotta clear visited since we didn't loop */
        memset(vctx->visited, 0, ctx->img_w * ctx->img_h * sizeof(uint8_t));
        memset(vctx->visited_m, 0, ctx->img_w * ctx->img_h * sizeof(uint8_t));
    }
}


void xy_sample(struct render_ctx *ctx, int x, int y, uint32_t n, uint32_t m, struct visited_ctx *vctx) {


    int o = y * ctx->img_w + x;


    for (uint32_t i = 0; i < n; i++) {

        uint32_t gcount = __sync_add_and_fetch(&(ctx->grid[o].count), 0);

        if (gcount < m) {

            __complex128 p = point_from_xy_rand(ctx, x, y);

            point_sample(ctx, p, vctx);
        } else {
            break;
        }
    }
}


void * image_sample_thread(void *targ) {

    struct thread_ctx *tctx = (struct thread_ctx *)targ;
    struct render_ctx *ctx = tctx->ctx;

    struct visited_ctx svctx;
    struct visited_ctx *vctx = &svctx;

    vctx->limit = ORD_LIMIT;

    vctx->visited = calloc(ctx->img_w * ctx->img_h, sizeof(uint8_t));
    assert(vctx->visited != NULL);
    vctx->visited_m = calloc(ctx->img_w * ctx->img_h, sizeof(uint8_t));
    assert(vctx->visited_m != NULL);

    vctx->vx = malloc(GROW_VISITED * sizeof(uint32_t));
    assert(vctx->vx != NULL);
    vctx->vy = malloc(GROW_VISITED * sizeof(uint32_t));
    assert(vctx->vy != NULL);

    vctx->vx_m = malloc(GROW_VISITED * sizeof(uint32_t));
    assert(vctx->vx_m != NULL);
    vctx->vy_m = malloc(GROW_VISITED * sizeof(uint32_t));
    assert(vctx->vy_m != NULL);

    vctx->vused = 0;
    vctx->vsize = GROW_VISITED;

    vctx->vused_m = 0;
    vctx->vsize_m = GROW_VISITED;


    fprintf(stderr, "== SAMPLING BORDER PIXELS ==\n");
    for (int y = tctx->tnum; y < ctx->img_h; y += tctx->num_threads) {

        fprintf(stderr, "Working on row %d of %d\n", y, ctx->img_h);

        for (int x = 0; x < ctx->img_w; x++) {

            if (xy_on_border(ctx, x, y) == 1) {
                xy_sample(ctx, x, y, 256, 128, vctx);
            }
        }
    }

    fprintf(stderr, "== SAMPLING DISC PIXELS ==\n");
    for (int y = tctx->tnum; y < ctx->img_h; y += tctx->num_threads) {
        fprintf(stderr, "Working on row %d of %d\n", y, ctx->img_h);

        for (int x = 0; x < ctx->img_w; x++) {

            xy_sample(ctx, x, y, 128, 64, vctx);
        }
    }


    free(vctx->visited);
    free(vctx->visited_m);
    free(vctx->vx);
    free(vctx->vy);
    free(vctx->vx_m);
    free(vctx->vy_m);

    return NULL;
}


void test_xy_point(struct render_ctx *ctx) {

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
                double log_avg_order = fabs((double)ctx->grid[o].scaled_log_order /
                                            ((double)ctx->grid[o].count * (double)LOG_SCALE));
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

    /* double max_order = 0; */
    /* double min_order = ORD_LIMIT; */
    /* for (int y = 0; y < ctx->img_h; y++) { */
    /*     for (int x = 0; x < ctx->img_w; x++) { */
    /*         int o = y * ctx->img_w + x; */
    /*         if (ctx->grid[o].count > 0) { */

    /*             double order; */
    /*             if (ctx->grid[o].ord_a >= ctx->grid[o].ord_b) { */
    /*                 order = (double)ctx->grid[o].ord_a / (double)ctx->grid[o].ord_b; */
    /*             } else { */
    /*                 order = (double)ctx->grid[o].ord_b / (double)ctx->grid[o].ord_a; */
    /*             } */

    /*             double avg_order = order / (double)ctx->grid[o].count; */
    /*             if (avg_order > max_order) { */
    /*                 max_order = avg_order; */
    /*             } */
    /*             if (avg_order < min_order) { */
    /*                 min_order = avg_order; */
    /*             } */
    /*         } */
    /*     } */
    /* } */


    /* color pixles scaled by max order */
    for (int y = 0; y < ctx->img_h; y++) {
        for (int x = 0; x < ctx->img_w; x++) {
            int o = y * ctx->img_w + x;
            if (ctx->grid[o].count > 0) {

                /* Figure out if we should blend with black */
                double bright_factor = 1.0;
                if (xy_on_border(ctx, x, y) == 1) {

                    int count_in_puzzle = 0;
                    for (int i = 0; i < 256; i++) {

                        /* Test point -> xy and xy -> point */
                        int nx, ny;
                        __complex128 tp;

                        nx = 0; ny = 0;

                        tp = point_from_xy_rand(ctx, x, y);

                        assert(point_to_xy(ctx, tp, &nx, &ny) == 0);

                        if ((nx != x) || (ny != y)) {
                            /* If this happens there is a bug with generating points
                             * in the desired x/y pixel */
                            fprintf(stderr, "Failed for point (%.5f, %.5f)\n", (double)__real__ tp, (double)__imag__ tp);
                        }

                        if (point_in_puzzle(ctx, tp) == 1) {
                            count_in_puzzle++;
                        }
                    }

                    bright_factor = (double)count_in_puzzle / 256.0;
                }


                int neg = 0;

                /* double order; */
                /* if (ctx->grid[o].ord_a >= ctx->grid[o].ord_b) { */
                /*     order = (double)ctx->grid[o].ord_a / (double)ctx->grid[o].ord_b; */
                /* } else { */
                /*     neg = 1; */
                /*     order = (double)ctx->grid[o].ord_b / (double)ctx->grid[o].ord_a; */
                /* } */

                /* double avg_order = order / (double)ctx->grid[o].count; */

                double log_avg_order = ((double)ctx->grid[o].scaled_log_order /
                                        ((double)ctx->grid[o].count * (double)LOG_SCALE));

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
                    image_data[o * 3 + 0] = (int)round(sin(v) * 255.0 * bright_factor);
                    image_data[o * 3 + 1] = (int)round((1.0 - sin(v * 2.0)) * 255.0 * bright_factor);
                    image_data[o * 3 + 2] = (int)round(cos(v) * 255.0 * bright_factor);
                } else {
                    image_data[o * 3 + 0] = (int)round((1.0 - sin(v)) * 255.0 * bright_factor);
                    image_data[o * 3 + 1] = (int)round(sin(v * 2.0) * 255.0 * bright_factor);
                    image_data[o * 3 + 2] = (int)round((1.0 - cos(v)) * 255.0 * bright_factor);
                }
            }
        }
    }

    write_png_file(name, ctx->img_w, ctx->img_h, image_data);

    free(image_data);
}


int main (void) {

    struct render_ctx *ctx = calloc(1, sizeof(struct render_ctx));
    ctx->n = 12;
    ctx->r = sqrt(2.0);
    ctx->r_sq = 2.0;
    ctx->epsilon = 1e-16Q;

    double goalw = 512;
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
    ctx->pradius = sqrt(pow(ctx->half_pwidth, 2.0) + pow(ctx->half_pheight, 2.0));
    ctx->grid = calloc(ctx->img_w * ctx->img_h, sizeof(struct samples));

    pthread_mutex_t grid_mutex = PTHREAD_MUTEX_INITIALIZER;
    ctx->grid_mutex = &(grid_mutex);

    __real__ ctx->rot_a = cosq(M_PIq * ((__float128)2 / (__float128)ctx->n));
    __imag__ ctx->rot_a = sinq(M_PIq * ((__float128)2 / (__float128)ctx->n));
    ctx->rot_b = conjq(ctx->rot_a);

    fprintf(stderr, "Image size %dx%d\n", ctx->img_w, ctx->img_h);

    fprintf(stderr, "Image box (%.5f, %.5f) to (%.5f, %.5f)\n",
            ctx->xmin, ctx->ymin,
            ctx->xmax, ctx->ymax);

    struct thread_ctx *tctxs = (struct thread_ctx *)calloc(NUM_THREADS, sizeof(struct thread_ctx));
    for (int i = 0; i < NUM_THREADS; i++) {
        tctxs[i].tnum = i;
        tctxs[i].num_threads = NUM_THREADS;
        tctxs[i].ctx = ctx;
    }

    int err;

    for (int i = 0; i < NUM_THREADS; i++) {
        err = pthread_create(&(tctxs[i].tid), NULL, image_sample_thread, &(tctxs[i]));
        if (err != 0) {
            perror("error creating stats thread");
            exit(255);
        }
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(tctxs[i].tid, NULL);
    }

    ctx_to_png(ctx, "/tmp/test.png");


    free(ctx->grid);
    free(ctx);

    return 0;
}

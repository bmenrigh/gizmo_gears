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


#define ORD_LIMIT (20 * 1000 * 1000)
#define NUM_THREADS 24
#define LOG_SCALE 0x1000000
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
    __float128 r;
    __float128 r_sq;
    int img_w, img_h;
    double xmin, xmax, ymin, ymax;
    double pwidth, pheight, half_pwidth, half_pheight, pradius;
    __complex128 rot_a, rot_b;
    __float128 epsilon;
    struct samples *grid;
    pthread_mutex_t *grid_mutex;
    int wedge_only;
    int box_only;
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


int xy_to_offset(struct render_ctx *ctx, int x, int y) {
    return y * ctx->img_w + x;
}


int point_to_xy(struct render_ctx *ctx, __complex128 p, int *x, int *y) {

    double px = (double)(__real__ p);
    double py = (double)(__imag__ p);

    if ((px < ctx->xmin) || (px > ctx->xmax) ||
        (py < ctx->ymin) || (py > ctx->ymax)) {
        return -1;
    }

    *x = (int)floor((px - ctx->xmin) / ctx->pwidth);
    *y = (int)floor((ctx->ymax - py) / ctx->pheight);

    return 0;
}


__complex128 point_from_xy(struct render_ctx *ctx, int x, int y) {

    __complex128 p;

    __real__ p = ctx->xmin + (ctx->pwidth * x);
    __imag__ p = ctx->ymax - (ctx->pheight * y);

    return p;
}


__complex128 point_rand_offset(struct render_ctx *ctx) {

    __complex128 p;
    __real__ p = (__float128)(((double)rand() / (double)0x80000000) * ctx->pwidth);
    __imag__ p = 0.0Q - ((__float128)(((double)rand() / (double)0x80000000) * ctx->pheight));

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

    if (fabsq(__real__ p - __real__ q) + fabsq(__imag__ p - __imag__ q)  < ctx->epsilon) {
        return 1;
    } else {
        return 0;
    }
}


int point_in_a(struct render_ctx *ctx, __complex128 p) {

    __complex128 np = p;
    __real__ np += 1.0Q;

    /* Check triangle inequality first */
    if (fabsq(__real__ np) + fabsq(__imag__ np) < ctx->r) {
        return 1;

        /* Else check squared pythagorean */
    } else if (((__real__ np) * (__real__ np)) +
                ((__imag__ np) * (__imag__ np)) < ctx->r_sq) {
        return 1;
    } else {
        return 0;
    }
}


int point_in_b(struct render_ctx *ctx, __complex128 p) {

    __complex128 np = p;
    __real__ np -= 1.0Q;

    /* Check triangle inequality first */
    if (fabsq(__real__ np) + fabsq(__imag__ np) < ctx->r) {
        return 1;

        /* Else check squared pythagorean */
    } else if (((__real__ np) * (__real__ np)) +
               ((__imag__ np) * (__imag__ np)) < ctx->r_sq) {
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
    __real__ np += 1.0Q;

    dist = (double)cabsq(np);
    if ((dist <= (double)ctx->r + (2.0 * ctx->pradius)) && (dist >= (double)ctx->r - (2.0 * ctx->pradius))) {
        return 1;
    }

    /* Check distace to B */
    np = p;
    __real__ np -= 1.0Q;

    dist = (double)cabsq(np);
    if ((dist <= (double)ctx->r + (2.0 * ctx->pradius)) && (dist >= (double)ctx->r - (2.0 * ctx->pradius))) {
        return 1;
    }

    return 0;
}


int point_in_puzzle(struct render_ctx *ctx, __complex128 p) {

    int in_a = point_in_a(ctx, p);
    int in_b = point_in_b(ctx, p);

    if (ctx->wedge_only == 1) {
        if ((in_a == 1) && (in_b == 1)) {
            return 1;
        }
    } else {
        if ((in_a == 1) || (in_b == 1)) {
            return 1;
        }
    }

    return 0;
}


int point_in_wedge(struct render_ctx *ctx, __complex128 p) {

    if ((point_in_a(ctx, p) == 1) && (point_in_b(ctx, p) == 1)) {
        return 1;
    } else {
        return 0;
    }
}


int point_in_box(struct render_ctx *ctx, __complex128 p) {

    double px = (double)(__real__ p);
    double py = (double)(__imag__ p);

    if ((px > ctx->xmin) && (px < ctx->xmax) &&
        (py > ctx->ymin) && (py < ctx->ymax)) {
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

    if ((x >= 0) && (y >= 0)) {
        uint32_t o = xy_to_offset(ctx, x, y);
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
    }

    if ((x_m >= 0) && (y_m >= 0)) {
        uint32_t o_m = xy_to_offset(ctx, x_m, y_m);
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

            /* Only track before a is done */
            x = -1;
            y = -1;
            x_m = -1;
            y_m = -1;
            if (((ctx->wedge_only == 0) && (ctx->box_only == 0)) ||
                ((ctx->wedge_only == 1) && (point_in_wedge(ctx, p) == 1)) ||
                ((ctx->box_only == 1) && (point_in_box(ctx, p) == 1))) {

                assert(point_to_xy(ctx, p, &x, &y) == 0);
            }

            __complex128 p_m;
            __real__ p_m = 0.0Q - __real__ p;
            __imag__ p_m = 0.0Q - __imag__ p;

            if (((ctx->wedge_only == 0) && (ctx->box_only == 0)) ||
                ((ctx->wedge_only == 1) && (point_in_wedge(ctx, p_m) == 1)) ||
                ((ctx->box_only == 1) && (point_in_box(ctx, p_m) == 1))) {

                assert(point_to_xy(ctx, p_m, &x_m, &y_m) == 0);
            }

            /* Now add */
            add_visited_xy(ctx, vctx, x, y, x_m, y_m);

            /* Actually try to do A turn */
            if (point_in_a(ctx, p) == 1) {
                p = turn_a(ctx, p);
                count_a++;
            }
        } else {
            /* Try B turn */
            if (point_in_b(ctx, p) == 1) {
                p = turn_b(ctx, p);
                count_b++;
            }
        }

        step ^= 1; /* toggle between a and b */
        count++;

        if (count > vctx->limit) {

            /* Try to salvage this point if it's extremely close to 0 */
            /* Within one loop around of each other */
            if (abs(count_a - count_b) <= ctx->n) {
                /*fprintf(stderr, "Salvaged point at limit\n");*/
                return 0.0;
            } else if (count > ORD_LIMIT) {
                /*fprintf(stderr, "count over limit assuming 0\n");*/
                return 0.0; /* Just assume they were equal as a speedup hack */
            } else {
                /*fprintf(stderr, "count over limit but not at max, returning NAN\n");*/
                return NAN; /* We didn't try long enough to be sure */
            }

        }

        if (count % 10000000 == 0) {
            fprintf(stderr, "Done %d turns\n", count);
        }

        /*
         * Note the condition (step != 0) here is critical for
         * correctness.  There are sometimes points that cycle around
         * back to themselves but they return to their original
         * position with a finaly A' turn so the next turn would be
         * B. This is different than how they started which was to
         * start with the next turn (first turn) being A'.
         *
         * This has the effect of causing two points to be measured
         * with two different orders depending on where they happen to
         * be sampled in their orbit.
         *
         * For example with N = 12; R = sqrt(2): The point
         * -0.199445936742669, 0.736391888282075 returns back to
         * itself with the opposite turn parity from how it started.
         * If the orbit is terminated at this step the order is 62
         * turns.  However if another point in the orbit is sampled:
         * -0.210794173119580, 0.770061921590612 this point will cycle
         * through the point -0.199445936742669, 0.736391888282075 but
         * has an order of 216 turns.
         *
         * So if all the points in the orbit of either of these get
         * labeled with their order then the pixel containing
         * -0.199445936742669, 0.736391888282075 can recieve two
         * different colors.
         *
         * By enforcing step != 0 to continue the loop points must
         * return back to their original spot with the same turn
         * parity and this causes the order to be measured the same no
         * matter what point in the orbit is sampled.
         *
         * It took me nearly 20 hours over three days to find this
         * bug.
         *
         */
    } while ((step != 0) || (count < ctx->n) || (point_equal_epsilon(ctx, op, p) != 1));

    /* delta order a/b stuff */
    if ((count_a > 0) && (count_b > 0)) {

        /* int excess = (((count_a % ctx->n) - (count_b % ctx->n)) + ctx->n) % ctx->n; */

        /* if (excess != 0) { */
        /*     fprintf(stderr, "point (%.15f, %.15f) with order %d with %d, %d\n", (double)__real__ op, (double)__imag__ op, count, count_a, count_b); */
        /* } */

        return (log((double)count_a) - log((double)count_b)); /* log(a/b) */
    } else {
        /*fprintf(stderr, "a or b was zero, returning NAN\n");*/
        return NAN;
    }
}


void point_sample(struct render_ctx *ctx, __complex128 p, struct visited_ctx *vctx) {

    if (point_in_puzzle(ctx, p) != 1) {
        return;
    }

    if ((ctx->wedge_only == 1) && (point_in_wedge(ctx, p) != 1)) {
        return;
    }

    if ((ctx->box_only == 1) && (point_in_box(ctx, p) != 1)) {
        return;
    }

    /* Don't sample from bottom, mirror and sample from top instead */
    /* if (__imag__ p < 0.0Q) { */
    /*     __real__ p = 0.0Q - __real__ p; */
    /*     __imag__ p = 0.0Q - __imag__ p; */
    /* } */

    /* Don't sample from wedge, rotate A by half and start there */
    /* if (point_in_wedge(ctx, p) == 1) { */
    /*     __real__ p = 0.0Q - (2.0Q + __real__ p); */
    /*     __imag__ p = 0.0Q - __imag__ p; */
    /* } */

    vctx->vused = 0;
    vctx->vused_m = 0;

    double ord = point_order(ctx, p, vctx);

    if (isnan(ord) == 0) {

        int64_t scaled_ord = (int64_t)round(ord * (double)LOG_SCALE);

        for (int i = 0; i < vctx->vused; i++) {
            int o = xy_to_offset(ctx, vctx->vx[i], vctx->vy[i]);

            __sync_add_and_fetch(&(ctx->grid[o].count), vctx->visited[o]);
            __sync_add_and_fetch(&(ctx->grid[o].scaled_log_order), scaled_ord * vctx->visited[o]);

            vctx->visited[o] = 0; /* clear this visit */
        }

        for (int i = 0; i < vctx->vused_m; i++) {
            int o_m = xy_to_offset(ctx, vctx->vx_m[i], vctx->vy_m[i]);

            __sync_add_and_fetch(&(ctx->grid[o_m].count), vctx->visited_m[o_m]);
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

    int o = xy_to_offset(ctx, x, y);

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
                xy_sample(ctx, x, y, 128, 128, vctx);
            }
        }
    }

    fprintf(stderr, "== SAMPLING DISC PIXELS ==\n");
    for (int y = tctx->tnum; y < ctx->img_h; y += tctx->num_threads) {

        fprintf(stderr, "Working on row %d of %d\n", y, ctx->img_h);

        for (int x = 0; x < ctx->img_w; x++) {
            xy_sample(ctx, x, y, 32, 4, vctx);
        }
    }

    /* fprintf(stderr, "== OVER SAMPLING LOW-ORDER PIXELS ==\n"); */
    /* vctx->limit = 100000; */
    /* for (int y = tctx->tnum; y < ctx->img_h; y += tctx->num_threads) { */
    /*     fprintf(stderr, "Working on row %d of %d\n", y, ctx->img_h); */

    /*     for (int x = 0; x < ctx->img_w; x++) { */

    /*         xy_sample(ctx, x, y, 128, 1, vctx); */
    /*     } */
    /* } */


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
            int o = xy_to_offset(ctx, x, y);
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

    double max_order = 11.0;
    double min_order = 1.0;
    /* double max_order = exp(log_max_order); */
    /* double min_order = exp(log_min_order); */
    fprintf(stderr, "log min: %f; min: %f\n", log_min_order, min_order);
    fprintf(stderr, "log max: %f; max: %f\n", log_max_order, max_order);


    /* color pixles scaled by max order */
    for (int y = 0; y < ctx->img_h; y++) {
        for (int x = 0; x < ctx->img_w; x++) {
            int o = xy_to_offset(ctx, x, y);
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

                double log_avg_order = ((double)ctx->grid[o].scaled_log_order /
                                        ((double)ctx->grid[o].count * (double)LOG_SCALE));

                if (log_avg_order < 0) {
                    neg = 1;
                    log_avg_order = fabs(log_avg_order);
                }

                double avg_order = exp(log_avg_order);
                double v = atan2(avg_order - min_order, 1.0) / atan2(max_order - min_order, 1.0);
                if (v > 1.0) {
                    v = 1.0;
                }

                /* delta colors */
                double vp2 = v * M_PI_2;

                if (neg == 0) {
                    image_data[o * 3 + 0] = (int)round(sin(vp2) * 255.0 * bright_factor);
                    image_data[o * 3 + 1] = (int)round((1.0 - sin(vp2 * 2.0)) * 255.0 * bright_factor);
                    image_data[o * 3 + 2] = (int)round(cos(vp2) * 255.0 * bright_factor);
                } else {
                    image_data[o * 3 + 0] = (int)round((1.0 - sin(vp2)) * 255.0 * bright_factor);
                    image_data[o * 3 + 1] = (int)round(sin(vp2 * 2.0) * 255.0 * bright_factor);
                    image_data[o * 3 + 2] = (int)round((1.0 - cos(vp2)) * 255.0 * bright_factor);
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
    ctx->r = sqrtq(2.0Q);
    ctx->r_sq = 2.0Q;
    ctx->epsilon = 1e-16Q;

    ctx->wedge_only = 1;
    ctx->box_only = 0;

    /* Render full puzzle */
    /* double goalw = 1024; */
    /* double scalef = goalw / ((2.0 * ctx->r) + 2.0); */
    /* ctx->img_w = (int)floor(((2.0 * ctx->r) + 2.0) * scalef); */
    /* ctx->img_h = (int)floor(2.0 * ctx->r * scalef); */
    /* ctx->xmin = -1.0 - ctx->r; */
    /* ctx->xmax = 1.0 + ctx->r; */
    /* ctx->ymin = -1.0 * ctx->r; */
    /* ctx->ymax = 1.0 * ctx->r; */

    /* Render wedge only */
    double goalh = 1024;
    double wedge_height = sqrt(ctx->r_sq - 1.0);
    double wedge_width = ctx->r - 1.0;
    ctx->img_h = (int)floor(goalh);
    ctx->img_w = (int)floor(goalh * (wedge_width / wedge_height));
    ctx->xmin = 0.0 - wedge_width;
    ctx->xmax = wedge_width;
    ctx->ymin = 0.0 - wedge_height;
    ctx->ymax = wedge_height;

    /* Render box only */
    /* double goalw = 2048; */
    /* ctx->xmin = -0.198298374860865 - 0.05; */
    /* ctx->xmax = -0.198298374860865 + 0.05; */
    /* ctx->ymin = 0.735207472673922 - 0.05; */
    /* ctx->ymax = 0.735207472673922 + 0.05; */
    /* ctx->img_w = (int)floor(goalw); */
    /* ctx->img_h = (int)floor(goalw / ((ctx->xmax - ctx->xmin) / (ctx->ymax - ctx->ymin))); */

    ctx->pwidth = (ctx->xmax - ctx->xmin) / (double)ctx->img_w;
    ctx->pheight = (ctx->ymax - ctx->ymin) / (double)ctx->img_h;
    ctx->half_pwidth = ctx->pwidth / 2.0;
    ctx->half_pheight = ctx->pheight / 2.0;
    ctx->pradius = sqrt(pow(ctx->half_pwidth, 2.0) + pow(ctx->half_pheight, 2.0));
    ctx->grid = calloc(ctx->img_w * ctx->img_h, sizeof(struct samples));

    pthread_mutex_t grid_mutex = PTHREAD_MUTEX_INITIALIZER;
    ctx->grid_mutex = &(grid_mutex);

    __real__ ctx->rot_a = cosq(M_PIq * (2.0Q / (__float128)ctx->n));
    __imag__ ctx->rot_a = sinq(M_PIq * (2.0Q / (__float128)ctx->n));
    ctx->rot_b = conjq(ctx->rot_a);

    test_xy_point(ctx);

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

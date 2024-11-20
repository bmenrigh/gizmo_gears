
#define  _GNU_SOURCE /* gets us M_PIl */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>
#include <quadmath.h>
#include <complex.h>
#include <assert.h>
#include <setjmp.h>
#include <pthread.h>

#include <getopt.h>

#include <png.h>


#define ORD_LIMIT    (500 * 1000 * 1000)
#define NUM_THREADS  24
#define LOG_SCALE    0x100000000
#define GROW_VISITED 1024
#define MAX_GEN_LEN  8
#define FNAME_LEN    1024 /* max output file path length */

/* #define naming conventions:
 * _T -> type
 * _F -> function
 * _L -> literal
 */

#define FPREC 64

#if FPREC == 128

#define COMPLEX_T __complex128
#define FLOAT_T __float128
#define FLOAT_L(N) (N##Q)
#define FABS_F(x) (fabsq(x))
#define CABS_F(x) (cabsq(x))
#define SQRT_F(x) (sqrtq(x))
#define FSIN_F(x) (sinq(x))
#define FCOS_F(x) (cosq(x))
#define STRTO_F(x) (strtoflt128(x, NULL))
#define PI_L      M_PIq
#define EPSILON_L FLOAT_L(1e-16)
#define MANTDIG_L FLT128_MANT_DIG

#elif FPREC == 80

#define COMPLEX_T _Complex long double
#define FLOAT_T long double
#define FLOAT_L(N) (N##L)
#define FABS_F(x) (fabsl(x))
#define CABS_F(x) (cabsl(x))
#define SQRT_F(x) (sqrtl(x))
#define FSIN_F(x) (sinl(x))
#define FCOS_F(x) (cosl(x))
#define STRTO_F(x) (strtold(x, NULL))
#define PI_L      M_PIl
#define EPSILON_L FLOAT_L(1e-10)
#define MANTDIG_L LDBL_MANT_DIG

#elif FPREC == 64

#define COMPLEX_T complex
#define FLOAT_T double
#define FLOAT_L(N) (N)
#define FABS_F(x) (fabs(x))
#define CABS_F(x) (cabs(x))
#define SQRT_F(x) (sqrt(x))
#define FSIN_F(x) (sin(x))
#define FCOS_F(x) (cos(x))
#define STRTO_F(x) (strtod(x, NULL))
#define PI_L      M_PI
#define EPSILON_L FLOAT_L(1e-10)
#define MANTDIG_L DBL_MANT_DIG

# else

#error "FPREC must defined and set to one of {128,80,64}"

#endif


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
};

struct render_ctx {
    char out[FNAME_LEN];
    int verbose;
    uint32_t n;
    FLOAT_T r;
    FLOAT_T r_sq;
    int img_w, img_h;
    double xmin, xmax, ymin, ymax;
    double pwidth, pheight, half_pwidth, half_pheight, pradius;
    COMPLEX_T rot[MAX_GEN_LEN]; /* alternates a, b, a, b, ... */
    int gen_len;
    FLOAT_T epsilon;
    struct samples *grid;
    pthread_mutex_t *grid_mutex;
    int wedge_only;
    int box_only;
    int sym180;
    int order_delta;
    uint64_t highest_order;
    double stretch_exp; /* order to val stretching */
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


COMPLEX_T turn_angle(struct render_ctx *ctx, int amt) {

    COMPLEX_T rot_ang;

    __real__ rot_ang = FCOS_F(PI_L * FLOAT_L(-2.0) * ((FLOAT_T)amt / (FLOAT_T)ctx->n));
    __imag__ rot_ang = FSIN_F(PI_L * FLOAT_L(-2.0) * ((FLOAT_T)amt / (FLOAT_T)ctx->n));

    return rot_ang;
}


int xy_to_offset(struct render_ctx *ctx, int x, int y) {
    return y * ctx->img_w + x;
}


int point_to_xy(struct render_ctx *ctx, COMPLEX_T p, int *x, int *y) {

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


COMPLEX_T point_from_xy(struct render_ctx *ctx, int x, int y) {

    COMPLEX_T p;

    __real__ p = ctx->xmin + (ctx->pwidth * x);
    __imag__ p = ctx->ymax - (ctx->pheight * y);

    return p;
}


COMPLEX_T point_rand_offset(struct render_ctx *ctx) {

    COMPLEX_T p;
    __real__ p = (FLOAT_T)(((double)rand() / (double)0x80000000) * ctx->pwidth);
    __imag__ p = FLOAT_L(0.0) - ((FLOAT_T)(((double)rand() / (double)0x80000000) * ctx->pheight));

    return p;
}


COMPLEX_T point_from_xy_rand(struct render_ctx *ctx, int x, int y) {

    return point_from_xy(ctx, x, y) + point_rand_offset(ctx);
}


int point_equal_double(COMPLEX_T p, COMPLEX_T q) {

    if (((double)(__real__ p) == (double)(__real__ q)) &&
        ((double)(__imag__ p) == (double)(__imag__ q))) {

        return 1;
    } else {
        return 0;
    }
}


int point_equal_epsilon(struct render_ctx *ctx, COMPLEX_T p, COMPLEX_T q) {

    if (FABS_F(__real__ p - __real__ q) + FABS_F(__imag__ p - __imag__ q)  < ctx->epsilon) {
        return 1;
    } else {
        return 0;
    }
}


int point_in_a(struct render_ctx *ctx, COMPLEX_T p) {

    COMPLEX_T np = p;
    __real__ np += FLOAT_L(1.0);

    /* Check triangle inequality first */
    if (FABS_F(__real__ np) + FABS_F(__imag__ np) < ctx->r) {
        return 1;

        /* Else check squared pythagorean */
    } else if (((__real__ np) * (__real__ np)) +
                ((__imag__ np) * (__imag__ np)) < ctx->r_sq) {
        return 1;
    } else {
        return 0;
    }
}


int point_in_n(struct render_ctx *ctx, COMPLEX_T p, int n) {

    COMPLEX_T np = p;

    if ((n & 1) == 0) {
        /* an A check */
        __real__ np += FLOAT_L(1.0);
    } else {
        /* a B check */
        __real__ np -= FLOAT_L(1.0);
    }

    /* Check triangle inequality first */
    if (FABS_F(__real__ np) + FABS_F(__imag__ np) < ctx->r) {
        return 1;

        /* Else check squared pythagorean */
    } else if (((__real__ np) * (__real__ np)) +
               ((__imag__ np) * (__imag__ np)) < ctx->r_sq) {
        return 1;
    } else {
        return 0;
    }
}


int point_in_b(struct render_ctx *ctx, COMPLEX_T p) {

    COMPLEX_T np = p;
    __real__ np -= FLOAT_L(1.0);

    /* Check triangle inequality first */
    if (FABS_F(__real__ np) + FABS_F(__imag__ np) < ctx->r) {
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

    COMPLEX_T p = point_from_xy(ctx, x, y);
    COMPLEX_T np;
    double dist;

    /* Check distace to A */
    np = p;
    __real__ np += FLOAT_L(1.0);

    dist = (double)CABS_F(np);
    if ((dist <= (double)ctx->r + (2.0 * ctx->pradius)) && (dist >= (double)ctx->r - (2.0 * ctx->pradius))) {
        return 1;
    }

    /* Check distace to B */
    np = p;
    __real__ np -= FLOAT_L(1.0);

    dist = (double)CABS_F(np);
    if ((dist <= (double)ctx->r + (2.0 * ctx->pradius)) && (dist >= (double)ctx->r - (2.0 * ctx->pradius))) {
        return 1;
    }

    return 0;
}


int point_in_puzzle(struct render_ctx *ctx, COMPLEX_T p) {

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


int point_in_wedge(struct render_ctx *ctx, COMPLEX_T p) {

    if ((point_in_a(ctx, p) == 1) && (point_in_b(ctx, p) == 1)) {
        return 1;
    } else {
        return 0;
    }
}


int point_in_box(struct render_ctx *ctx, COMPLEX_T p) {

    double px = (double)(__real__ p);
    double py = (double)(__imag__ p);

    if ((px > ctx->xmin) && (px < ctx->xmax) &&
        (py > ctx->ymin) && (py < ctx->ymax)) {
        return 1;
    } else {
        return 0;
    }
}


COMPLEX_T turn_a(struct render_ctx *ctx, COMPLEX_T p) {

    COMPLEX_T np = p;
    __real__ np += FLOAT_L(1.0);

    np *= ctx->rot[0];

    __real__ np -= FLOAT_L(1.0);

    return np;
}


COMPLEX_T turn_b(struct render_ctx *ctx, COMPLEX_T p) {

    COMPLEX_T np = p;
    __real__ np -= FLOAT_L(1.0);

    np *= ctx->rot[1];

    __real__ np += FLOAT_L(1.0);

    return np;
}


COMPLEX_T turn_n(struct render_ctx *ctx, COMPLEX_T p, int n) {

    COMPLEX_T np = p;
    if ((n & 1) == 0) {
        /* an A turn */

        __real__ np += FLOAT_L(1.0);

        np *= ctx->rot[n];

        __real__ np -= FLOAT_L(1.0);
    } else {
        /* a B turn */

        __real__ np -= FLOAT_L(1.0);

        np *= ctx->rot[n];

        __real__ np += FLOAT_L(1.0);
    }

    return np;
}


double convolve_get_xy_val(double *source, int w, int h, double edge_val, int x, int y) {

    if ((x < 0) || (x >= w) ||
        (y < 0) || (y >= h)) {

        return edge_val;
    }

    return source[y * w + x];
}


void convolve_sobel(double *source, double *mag, double *angle, int w, int h, double edge_val) {

    double mh[3][3] = {{1.0, 0.0, -1.0}, {2.0, 0.0, -2.0}, {1.0, 0.0, -1.0}};
    double mv[3][3] = {{1.0, 2.0, 1.0}, {0.0, 0.0, 0.0}, {-1.0, -2.0, -1.0}};

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {

            double sob_h = 0.0;
            double sob_v = 0.0;

            for (int yo = -1; yo <= 1; yo++) {
                for (int xo = -1; xo <= 1; xo++) {
                    sob_h += mh[yo + 1][xo + 1] * convolve_get_xy_val(source, w, h, edge_val, x + xo, y + yo);
                    sob_v += mv[yo + 1][xo + 1] * convolve_get_xy_val(source, w, h, edge_val, x + xo, y + yo);
                }
            }

            /* magnitude uses pythagorean theorem */
            mag[y * w + x] = sqrt((sob_h * sob_h) + (sob_v * sob_v));

            /* angle is atan of rise / run */
            angle[y * w + x] = atan2(sob_v, sob_h);
        }
    }
}


double delta_log_order_to_val(double log_order, double min_order, double max_order) {

    /* Only let exp() operate on positive logs */
    int neg = 0;
    if (log_order < 0) {
        neg = 1;
        log_order = fabs(log_order);
    }

    /* Undo log */
    double order = exp(log_order);

    /* Saturate order in case of minor roundoff issues */
    if (order > max_order) {
        order = max_order;
    }
    if (order < min_order) {
        order = min_order;
    }

    /* atan -> val mapping */
    double v = atan2(order - min_order, 1.0) / atan2(max_order - min_order, 1.0);

    /* Saturate (shouldn't happen if atan2 behaves itself) */
    if (v > 1.0) {
        v = 1.0;
    }

    /* Flip val to [0, -1] if the log was negative */
    if (neg == 1) {
        v = 0.0 - v;
    }

    return v;
}


double log_order_to_val(double log_order, double min_order, double max_order, double stretch_exp) {

    double log_min_order = log(min_order);
    double log_max_order = log(max_order);

    double offset = 1.0; /* required so log(x < 1) doesn't turn negative */

    double v = ((log_order - log_min_order) + offset) / ((log_max_order - log_min_order) + offset);

    /* This stretches out the color space where more colors are allocated
     * to higher orders. An exponent of 1 is no stretch. 2.0 is a minor stretch.
     */
    return 1.0 - pow((1.0 - v), stretch_exp);
}


double log_order_atan_to_val(double log_order, double log_min_order, double log_max_order, double stretch_exp) {

    /* Saturate order in case of minor roundoff issues */
    if (log_order > log_max_order) {
        log_order = log_max_order;
    }
    if (log_order < log_min_order) {
        log_order = log_min_order;
    }

    /* atan -> val mapping */
    double v = atan2(log_order - log_min_order, 1.0) / atan2(log_max_order - log_min_order, 1.0);

    /* Saturate (shouldn't happen if atan2 behaves itself) */
    if (v > 1.0) {
        v = 1.0;
    }

    return 1.0 - pow((1.0 - v), stretch_exp);
}


double log_order_scale_cycle_to_val(double log_order) {

    double neg = 1;
    if (log_order < 0) {
        neg = -1;
        log_order = fabs(log_order);
    }

    double v = fmod(log_order, log(10.0)); /* fractional part of log base 10 */

    return v * neg;
}


void val_to_rgb(double val, uint8_t *R, uint8_t *G, uint8_t *B, double brightness) {

    assert((brightness >= 0.0) && (brightness <= 1.0));

    /* delta colors */
    double vp2 = val * M_PI_2;

    if (vp2 >= 0) {
        *R = (uint8_t)round(sin(vp2) * 255.0 * brightness);
        *G = (uint8_t)round((1.0 - sin(vp2 * 2.0)) * 255.0 * brightness);
        *B = (uint8_t)round(cos(vp2) * 255.0 * brightness);
    } else {
        vp2 = 0.0 - vp2; /* flip to positive */

        *R = (int)round((1.0 - sin(vp2)) * 255.0 * brightness);
        *G = (int)round(sin(vp2 * 2.0) * 255.0 * brightness);
        *B = (int)round((1.0 - cos(vp2)) * 255.0 * brightness);
    }
}


double hsv_to_rgb_f(int n, double h, double s, double v) {
    /* https://en.wikipedia.org/wiki/HSL_and_HSV */

    double k = fmod((double)n + h * 6.0, 6.0);

    double k_adj = fmax(0.0, fmin(fmin(k, 4 - k), 1.0));

    return v - (v * s * k_adj);
}


void hsv_to_rgb(double h, double s, double v, uint8_t *R, uint8_t *G, uint8_t *B) {

    *R = (uint8_t)round(hsv_to_rgb_f(5.0, h, s, v) * 255.0);
    *G = (uint8_t)round(hsv_to_rgb_f(3.0, h, s, v) * 255.0);
    *B = (uint8_t)round(hsv_to_rgb_f(1.0, h, s, v) * 255.0);
}


void val_to_rgb2(double val, uint8_t *R, uint8_t *G, uint8_t *B, double brightness) {

    assert((brightness >= 0.0) && (brightness <= 1.0));

    /* delta colors */
    double vp2 = val * M_PI_2;

    double h, s, v;

    if (vp2 >= 0) {
        h = 0;
    } else {
        vp2 = 0.0 - vp2; /* flip to positive */

        h = 0.5;
    }

    s = 0.25 + 0.5 * sin(2.0 * vp2);

    v = 0.75 + 0.25 * sin(vp2);
    v *= brightness;

    hsv_to_rgb(h, s, v, R, G, B);
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


double point_order(struct render_ctx *ctx, COMPLEX_T p, struct visited_ctx *vctx) {

    if (point_in_puzzle(ctx, p) != 1) {
        return 0;
    }

    int x, y, x_m, y_m;
    COMPLEX_T op = p; /* Original p */
    uint8_t step = 0;
    uint32_t count = 0;
    int32_t count_a = 0;
    int32_t count_b = 0;
    int stuck = 0;
    do {

        if (step == 0) {

            /* Do stuck checking */
            if (stuck == 1) {
                fprintf(stderr, "point (%.15f, %.15f) stuck at count %d\n", (double)__real__ op, (double)__imag__ op, count);

                return NAN;
            }
            stuck = 1;

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

            if (ctx->sym180 == 1) {
                COMPLEX_T p_m;
                __real__ p_m = FLOAT_L(0.0) - __real__ p;
                __imag__ p_m = FLOAT_L(0.0) - __imag__ p;

                if (((ctx->wedge_only == 0) && (ctx->box_only == 0)) ||
                    ((ctx->wedge_only == 1) && (point_in_wedge(ctx, p_m) == 1)) ||
                    ((ctx->box_only == 1) && (point_in_box(ctx, p_m) == 1))) {

                    assert(point_to_xy(ctx, p_m, &x_m, &y_m) == 0);
                }
            }

            /* Now add */
            add_visited_xy(ctx, vctx, x, y, x_m, y_m);
        }

        /* Try to do turn */
        if (point_in_n(ctx, p, step) == 1) {
            count++;
            stuck = 0;

            p = turn_n(ctx, p, step);

            if ((step & 1) == 0) {
                count_a++;
            } else {
                count_b++;
            }
        }

        /*step ^= 1;*/ /* toggle between a and b */
        step = (step + 1) % ctx->gen_len;

        if ((count > 0) && (count % 50000000 == 0)) {
            fprintf(stderr, "Done %d turns\n", count);
        }

        /*
         * Note the condition (step != 0) here is critical for
         * correctness.  There are sometimes points that cycle around
         * back to themselves, but they return to their original
         * position partway through the application of the generator.
         * For A' B that means they end with a final A' turn (so the
         * next turn would be B). This is different than how they
         * started which was to start with the next turn (first turn)
         * being A'.
         *
         * This has the effect of causing some points to be measured
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
         * labeled with their order, then the pixel containing
         * -0.199445936742669, 0.736391888282075 can receive two
         * different order values.
         *
         * By enforcing step != 0 to continue, the loop points must
         * return back to their original spot with the same turn
         * parity, and this causes the order to be measured the same
         * no matter what point in the orbit is sampled.
         *
         * It took me nearly 20 hours over three days to find this
         * bug.
         *
         */
    } while ((count < vctx->limit) && ((step != 0) || (count < ctx->n) || (point_equal_epsilon(ctx, op, p) != 1)));

    uint64_t order_sum = count_a + count_b;

    if (count < vctx->limit) {
        uint64_t max_order = __atomic_load_n(&(ctx->highest_order), __ATOMIC_RELAXED);

        if (order_sum > max_order) {

            fprintf(stderr, "New max order of %lu found for point (%.15f, %.15f)\n", order_sum, (double)__real__ op, (double)__imag__ op);
        }

        /* Now write new order back into highest_seen atomically */
        while ((order_sum > max_order) && (__atomic_compare_exchange_n(&(ctx->highest_order), &max_order, order_sum, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED) == true));
    }


    /* Figure out how to calculate order based on absolute or delta measure */
    if (ctx->order_delta == 0) {
        /* absolute order a + b stuff */

        if (count >= vctx->limit) {

            fprintf(stderr, "point (%.15f, %.15f) hit limit with order %lu\n", (double)__real__ op, (double)__imag__ op, order_sum);

            /*return NAN;*/
        }
        if (order_sum == 0) {
            return NAN;
        }

        return (log((double)(order_sum))); /* log(a + b) */

    } else {
        /* delta order a/b stuff */

        if (count >= vctx->limit) {

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

        if ((count_a > 0) && (count_b > 0)) {

            /* int excess = (((count_a % ctx->n) - (count_b % ctx->n)) + ctx->n) % ctx->n; */

            /* if (excess != 0) { */
            /*     fprintf(stderr, "point (%.15f, %.15f) with order %d with %d, %d\n", (double)__real__ op, (double)__imag__ op, count, count_a, count_b); */
            /* } */

            return (log((double)count_a) - log((double)count_b)); /* log(a/b) */
        } else {

            /* Check that if this is one of the the points in the
             * center of a disk that can only be turned by that disk
             * and aren't reached by the other disk.
             */
            if (count_a == (int32_t)ctx->n) {
                return log((double)count_a);
            } else if (count_b == (int32_t)ctx->n) {
                return 0.0 - log((double)count_b);
            } else {
                /*fprintf(stderr, "a or b was zero, returning NAN\n");*/
                return NAN;
            }
        }
    }
}


void point_sample(struct render_ctx *ctx, COMPLEX_T p, struct visited_ctx *vctx) {

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

        if (ctx->sym180 == 1) {
            for (int i = 0; i < vctx->vused_m; i++) {
                int o_m = xy_to_offset(ctx, vctx->vx_m[i], vctx->vy_m[i]);

                __sync_add_and_fetch(&(ctx->grid[o_m].count), vctx->visited_m[o_m]);

                if (ctx->order_delta == 0) {
                    __sync_add_and_fetch(&(ctx->grid[o_m].scaled_log_order), scaled_ord * vctx->visited_m[o_m]);
                } else {
                    __sync_sub_and_fetch(&(ctx->grid[o_m].scaled_log_order), scaled_ord * vctx->visited_m[o_m]);
                }

                vctx->visited_m[o_m] = 0; /* clear this visit */
            }
        }
    } else {
        /* Gotta clear visited since we didn't loop */
        memset(vctx->visited, 0, ctx->img_w * ctx->img_h * sizeof(uint8_t));

        if (ctx->sym180 == 1) {
            memset(vctx->visited_m, 0, ctx->img_w * ctx->img_h * sizeof(uint8_t));
        }
    }
}


void xy_sample(struct render_ctx *ctx, int x, int y, uint32_t n, uint32_t m, struct visited_ctx *vctx) {

    int o = xy_to_offset(ctx, x, y);

    for (uint32_t i = 0; i < n; i++) {

        uint32_t gcount = __sync_add_and_fetch(&(ctx->grid[o].count), 0);

        if (gcount < m) {

            COMPLEX_T p = point_from_xy_rand(ctx, x, y);

            point_sample(ctx, p, vctx);
        } else {
            break;
        }
    }
}


void image_aa_sobel(struct render_ctx *ctx) {

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


    fprintf(stderr, "== OVER SAMPLING SOBEL EDGE PIXELS ==\n");

    double max_order;
    double min_order;

    if (ctx->order_delta == 0) {
        max_order = vctx->limit;
        min_order = (double)ctx->n;
    } else {
        max_order = (double)ctx->n;
        min_order = 1.0;
    }

    double *vgrid = calloc(ctx->img_h * ctx->img_w, sizeof(double));
    double *vsobel_mag = calloc(ctx->img_h * ctx->img_w, sizeof(double));
    double *vsobel_ang = calloc(ctx->img_h * ctx->img_w, sizeof(double));

    /* Fill in the vgrid with vals for each pixel */
    for (int y = 0; y < ctx->img_h; y++) {
        for (int x = 0; x < ctx->img_w; x++) {

            int o = xy_to_offset(ctx, x, y);

            if (ctx->grid[o].count > 0) {

                double log_avg_order = ((double)ctx->grid[o].scaled_log_order /
                                        ((double)ctx->grid[o].count * (double)LOG_SCALE));

                double v;
                if (ctx->order_delta == 0) {
                    v = log_order_to_val(log_avg_order, min_order, max_order, ctx->stretch_exp);
                    /*v = log_order_atan_to_val(log_avg_order, min_order, max_order, ctx->stretch_exp);*/
                    /*v = log_order_scale_cycle_to_val(log_avg_order);*/
                } else {
                    v = delta_log_order_to_val(log_avg_order, min_order, max_order);
                }

                vgrid[o] = v;
            }
        }
    }

    /* Run sobel filter on vgrid */
    convolve_sobel(vgrid, vsobel_mag, vsobel_ang, ctx->img_w, ctx->img_h, 0.0);

    /* normalize magnitudes back into [-1, 1] */
    /* first find max maginutude */
    double max_mag = 1.0;
    for (int y = 0; y < ctx->img_h; y++) {
        for (int x = 0; x < ctx->img_w; x++) {

            int o = xy_to_offset(ctx, x, y);

            if (fabs(vsobel_mag[o]) > max_mag) {
                max_mag = fabs(vsobel_mag[o]);
            }
        }
    }

    /* now scale down by max */
    for (int y = 0; y < ctx->img_h; y++) {
        for (int x = 0; x < ctx->img_w; x++) {

            int o = xy_to_offset(ctx, x, y);

            vsobel_mag[o] /= max_mag;
        }
    }

    /* Use the sobel filter magnitude to select pixels for
     * sampling to make sure edges have enough samples
     */
    for (int y = 0; y < ctx->img_h; y++) {

        if (y % 25 == 0) {
            fprintf(stderr, "SOBEL: working on row %d of %d\n", y, ctx->img_h);
        }

        for (int x = 0; x < ctx->img_w; x++) {

            int o = xy_to_offset(ctx, x, y);
            int amt = (int)(1024.0 * fabs(vsobel_mag[o]));

            xy_sample(ctx, x, y, 1024, amt, vctx);
        }
    }


    free(vgrid);
    free(vsobel_mag);
    free(vsobel_ang);

    free(vctx->visited);
    free(vctx->visited_m);
    free(vctx->vx);
    free(vctx->vy);
    free(vctx->vx_m);
    free(vctx->vy_m);

    return;
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

    /* COMPLEX_T p; */
    /* __real__ p = 0; */
    /* __imag__ p = 0; */

    /* point_sample(ctx, p, vctx); */

    /* for (int i = 0; i < 100; i++) { */
    /*     xy_sample(ctx, (int)(ctx->img_w / 2), (int)(ctx->img_h / 2), 100, 100, vctx); */
    /* } */

    if (tctx->tnum == 0) {
        fprintf(stderr, "== SAMPLING BORDER PIXELS ==\n");
    }
    for (int y = tctx->tnum; y < ctx->img_h; y += tctx->num_threads) {

        if (tctx->tnum == 0) {
            fprintf(stderr, "BORDER: working on row %d of %d\n", y, ctx->img_h);
        }

        for (int x = 0; x < ctx->img_w; x++) {

            if (xy_on_border(ctx, x, y) == 1) {
                xy_sample(ctx, x, y, 256, 256, vctx);
            }
        }
    }

    if (tctx->tnum == 0) {
        fprintf(stderr, "== SAMPLING DISC PIXELS ==\n");
    }
    for (int y = tctx->tnum; y < ctx->img_h; y += tctx->num_threads) {

        if (tctx->tnum == 0) {
            fprintf(stderr, "DISC: working on row %d of %d\n", y, ctx->img_h);
        }

        for (int x = 0; x < ctx->img_w; x++) {
            xy_sample(ctx, x, y, 32, 4, vctx);
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
    COMPLEX_T tp;
    for (int i = 0; i < 1000; i++) {
        nx = 0; ny = 0;
        tp = point_from_xy_rand(ctx, tx, ty);
        assert(point_to_xy(ctx, tp, &nx, &ny) == 0);

        if ((nx != tx) || (ny != ty)) {
            fprintf(stderr, "Failed for point (%.5f, %.5f)\n", (double)__real__ tp, (double)__imag__ tp);
        }
    }

}


void ctx_to_png(struct render_ctx *ctx) {

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

    double max_order;
    double min_order;

    if (ctx->order_delta == 0) {
        max_order = exp(log_max_order);
        min_order = exp(log_min_order);
    } else {
        max_order = (double)ctx->n;
        min_order = 1.0;
    }

    fprintf(stderr, "log min: %f; min: %f\n", log_min_order, min_order);
    fprintf(stderr, "log max: %f; max: %f\n", log_max_order, max_order);

    double *vgrid = calloc(ctx->img_h * ctx->img_w, sizeof(double));
    double *vsobel_mag = calloc(ctx->img_h * ctx->img_w, sizeof(double));
    double *vsobel_ang = calloc(ctx->img_h * ctx->img_w, sizeof(double));

    /* Fill in the vgrid with vals for each pixel */
    for (int y = 0; y < ctx->img_h; y++) {
        for (int x = 0; x < ctx->img_w; x++) {

            int o = xy_to_offset(ctx, x, y);

            if (ctx->grid[o].count > 0) {

                double log_avg_order = ((double)ctx->grid[o].scaled_log_order /
                                        ((double)ctx->grid[o].count * (double)LOG_SCALE));

                double v;
                if (ctx->order_delta == 0) {
                    v = log_order_to_val(log_avg_order, min_order, max_order, ctx->stretch_exp);
                    /*v = log_order_atan_to_val(log_avg_order, log_min_order, log_max_order, ctx->stretch_exp);*/
                    /*v = log_order_scale_cycle_to_val(log_avg_order);*/
                } else {
                    v = delta_log_order_to_val(log_avg_order, min_order, max_order);
                }

                vgrid[o] = v;
            }
        }
    }

    /* Run sobel filter on vgrid */
    convolve_sobel(vgrid, vsobel_mag, vsobel_ang, ctx->img_w, ctx->img_h, 0.0);

    /* normalize magnitudes back into [-1, 1] */
    /* first find max maginutude */
    double max_mag = 1.0;
    for (int y = 0; y < ctx->img_h; y++) {
        for (int x = 0; x < ctx->img_w; x++) {

            int o = xy_to_offset(ctx, x, y);

            if (fabs(vsobel_mag[o]) > max_mag) {
                max_mag = fabs(vsobel_mag[o]);
            }
        }
    }

    /* now scale down by max */
    for (int y = 0; y < ctx->img_h; y++) {
        for (int x = 0; x < ctx->img_w; x++) {

            int o = xy_to_offset(ctx, x, y);

            vsobel_mag[o] /= max_mag;
        }
    }


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
                        COMPLEX_T tp;

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


                uint8_t *R = &(image_data[o * 3 + 0]);
                uint8_t *G = &(image_data[o * 3 + 1]);
                uint8_t *B = &(image_data[o * 3 + 2]);

                val_to_rgb(vgrid[o], R, G, B, bright_factor);
                /*val_to_rgb2(vgrid[o], R, G, B, bright_factor);*/
            }
        }
    }

    write_png_file(ctx->out, ctx->img_w, ctx->img_h, image_data);

    free(vgrid);
    free(vsobel_mag);
    free(vsobel_ang);

    free(image_data);
}


int main (int argc, char **argv) {

    struct render_ctx *ctx = calloc(1, sizeof(struct render_ctx));

    strncpy(ctx->out, "/tmp/test.png", FNAME_LEN - 1);

    ctx->highest_order = 0;
    ctx->epsilon = EPSILON_L;

    pthread_mutex_t grid_mutex = PTHREAD_MUTEX_INITIALIZER;
    ctx->grid_mutex = &(grid_mutex);

    ctx->wedge_only = 0;
    ctx->box_only = 0;
    ctx->order_delta = 0;

    ctx->n = 12;
    ctx->r = SQRT_F(2.0);
    ctx->r_sq = 2.0;
    ctx->sym180 = 1;

    ctx->stretch_exp = 1.0;

    int got_r = 0;
    int got_size = 0;
    int goalh = 512; /* default sets height */
    int goalw = 0; /* default sets height not width */

    static struct option long_options[] =
        {
            {"output",     required_argument,  0, 'o'},

            {"width",     required_argument,   0, 'x'},
            {"height",     required_argument,  0, 'y'},

            {"radius",     required_argument,  0, 'r'},
            {"n",          required_argument,  0, 'n'},
            {"radius-sq",  required_argument,  0,   0},

            {"color-delta", no_argument,       0,   0},
            {"stretch-exp", required_argument, 0,   0},
            {"wedge-only",  no_argument,       0,   0},

            {"verbose",     no_argument,       0, 'v'},

            /* Terminate list */
            {0, 0, 0, 0}
        };

    while (1) {

        /* getopt_long stores the option index here. */
        int option_index = 0;

        int c = getopt_long(argc, argv, "o:r:n:v:",
                            long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1) {
            break;
        }

        int ret;

        switch (c) {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0) {
                break;
            }

            if (strcmp("radius-sq", long_options[option_index].name) == 0) {
                ctx->r_sq = STRTO_F(optarg);
                ctx->r= SQRT_F(ctx->r_sq);
                got_r++;

                if (ctx->r_sq <= FLOAT_L(0.0)) {
                    fprintf(stderr, "Invalid radius-sq number found: \"%s\". Aborting.\n", optarg);
                    exit(-1);
                }

                break;
            }

            if (strcmp("color-delta", long_options[option_index].name) == 0) {
                ctx->order_delta = 1;
                break;
            }

            if (strcmp("stretch-exp", long_options[option_index].name) == 0) {
                ctx->stretch_exp = strtod(optarg, NULL);

                if (ctx->stretch_exp <= 0.0) {
                    fprintf(stderr, "Invalid stretch-exp number found: \"%s\". Aborting.\n", optarg);
                    exit(-1);
                }

                break;
            }


            if (strcmp("wedge-only", long_options[option_index].name) == 0) {
                ctx->wedge_only = 1;
                break;
            }

            if (strcmp("verbose", long_options[option_index].name) == 0) {
                ctx->verbose = 1;
                break;
            }

            abort(); /* should be unreachable */

        case 'x':
                ret = sscanf(optarg, "%u", &goalw);
                assert(ret == 1);
                goalh = 0;

                got_size++;

                if (goalw < 128) {
                    fprintf(stderr, "Width must be at least 128\n");
                    exit(-1);
                }

                break;

        case 'y':
                ret = sscanf(optarg, "%u", &goalh);
                assert(ret == 1);
                goalw = 0;

                got_size++;

                if (goalh < 128) {
                    fprintf(stderr, "Height must be at least 128\n");
                    exit(-1);
                }

                break;

        case 'o':
            strncpy(ctx->out, optarg, FNAME_LEN - 1);
            break;

        case 'r':
            ctx->r = STRTO_F(optarg);
            ctx->r_sq = ctx->r * ctx->r;
            got_r++;

            if (ctx->r <= FLOAT_L(0.0)) {
                fprintf(stderr, "Invalid radius number found: \"%s\". Aborting.\n", optarg);
                exit(-1);
            }

            break;

        case 'n':
            ret = sscanf(optarg, "%u", &(ctx->n));
            assert(ret == 1);
            break;

        case 'v':
            ctx->verbose = 1;
            break;

        default:
            fprintf(stderr, "Unrecognized option \"%c\". Aborting.\n", c);
            exit(-1);
        }
    }

    if (got_r > 1) {
        fprintf(stderr, "You can only set the radius once (-r and --radius-sq are mutually exclusive)\n");
        exit(-1);
    }

    if (got_size > 1) {
        fprintf(stderr, "You can only set the size once (-x and -y are mutually exclusive)\n");
        exit(-1);
    }


    /* Set image size based on specified goalh or goalw */
    double aspect_r;
    if (ctx->wedge_only == 0) {
        double discs_height = 2.0 * (double)ctx->r;
        double discs_width = (2.0 * (double)ctx->r) + 2.0;

        aspect_r = discs_width / discs_height;

        ctx->xmin = -1.0 - ctx->r;
        ctx->xmax = 1.0 + ctx->r;
        ctx->ymin = -1.0 * ctx->r;
        ctx->ymax = 1.0 * ctx->r;
    } else {
        double wedge_half_height = sqrt((double)ctx->r_sq - 1.0);
        double wedge_half_width = (double)ctx->r - 1.0;

        aspect_r = wedge_half_width / wedge_half_height;

        ctx->xmin = 0.0 - wedge_half_width;
        ctx->xmax = wedge_half_width;
        ctx->ymin = 0.0 - wedge_half_height;
        ctx->ymax = wedge_half_height;
    }

    if (goalh != 0) {
        goalw = (int)round((double)goalh * aspect_r);
    } else {
        goalh = (int)round((double)goalw / aspect_r);
    }

    ctx->img_w = goalw;
    ctx->img_h = goalh;


    /* ctx->n = 10; */
    /* ctx->r = SQRT_F((FLOAT_L(1.0) / FLOAT_L(2.0)) * (FLOAT_L(7.0) - SQRT_F(5.0))); */
    /* ctx->r_sq = (FLOAT_L(1.0) / FLOAT_L(2.0)) * (FLOAT_L(7.0) - SQRT_F(5.0)); */
    /* ctx->sym180 = 1; */

    /* ctx->n = 7; */
    /* ctx->r = FLOAT_L(1.62357492692335958); /\* tom rokicki approx *\/ */
    /* ctx->r = FLOAT_L(1.62318927293966937); /\* Eric's construction *\/ */
    /* ctx->r_sq = ctx->r * ctx->r; */
    /* ctx->sym180 = 0; */

    /* ctx->n = 12; */
    /* ctx->r = sqrtq(FLOAT_L(2.0)); */
    /* ctx->r_sq = FLOAT_L(2.0); */
    /* ctx->sym180 = 1; */

    /* n=12 critical radius */
    /* https://twistypuzzles.com/forum/viewtopic.php?t=25752&hilit=gizmo+gears+critical+radius&start=600 */
    /* benpuzzles: R=sqrt(40-22*sqrt(3)) or 1.376547214... */
    /* ctx->n = 12; */
    /* ctx->r = sqrtq(40.0Q - 22.0Q * sqrtq(3.0Q)); */
    /* ctx->r_sq = 40.0Q - 22.0Q * sqrtq(3.0Q); */

    /* N=5 critical radius */
    /* ctx->n = 5; */
    /* ctx->r = sqrtq((7.0Q + sqrtq(5.0Q)) / 2.0Q); */
    /* ctx->r_sq = (7.0Q + sqrtq(5.0Q)) / 2.0Q; */



    /* ctx->n = 14; */
    /* ctx->r_sq = 1.6180339887Q; */
    /* ctx->r = sqrtq(ctx->r_sq); */


    /* Render full puzzle */
    /* double goalw = 2048; */
    /* double scalef = goalw / ((2.0 * ctx->r) + 2.0); */
    /* ctx->img_w = (int)floor(((2.0 * ctx->r) + 2.0) * scalef); */
    /* ctx->img_h = (int)floor(2.0 * ctx->r * scalef); */
    /* ctx->xmin = -1.0 - ctx->r; */
    /* ctx->xmax = 1.0 + ctx->r; */
    /* ctx->ymin = -1.0 * ctx->r; */
    /* ctx->ymax = 1.0 * ctx->r; */

    /* Render wedge only */
    /* double goalh = 4096; */
    /* double wedge_height = sqrt((double)ctx->r_sq - 1.0); */
    /* double wedge_width = (double)ctx->r - 1.0; */
    /* ctx->img_h = (int)floor(goalh); */
    /* ctx->img_w = (int)floor(goalh * (wedge_width / wedge_height)); */
    /* ctx->xmin = 0.0 - wedge_width; */
    /* ctx->xmax = wedge_width; */
    /* ctx->ymin = 0.0 - wedge_height; */
    /* ctx->ymax = wedge_height; */

    /* Render box only */
    /* double goalw = 4096; */
    /* ctx->xmin = -0.1; */
    /* ctx->xmax = 0.0; */
    /* ctx->ymin = -0.1; */
    /* ctx->ymax = 0.0; */
    /* ctx->img_w = (int)floor(goalw); */
    /* ctx->img_h = (int)floor(goalw / ((ctx->xmax - ctx->xmin) / (ctx->ymax - ctx->ymin))); */


    ctx->pwidth = (ctx->xmax - ctx->xmin) / (double)ctx->img_w;
    ctx->pheight = (ctx->ymax - ctx->ymin) / (double)ctx->img_h;
    ctx->half_pwidth = ctx->pwidth / 2.0;
    ctx->half_pheight = ctx->pheight / 2.0;
    ctx->pradius = sqrt(pow(ctx->half_pwidth, 2.0) + pow(ctx->half_pheight, 2.0));
    ctx->grid = calloc(ctx->img_w * ctx->img_h, sizeof(struct samples));


    /* A' B generator */
    ctx->gen_len = 2;
    ctx->rot[0] = turn_angle(ctx, -1);
    ctx->rot[1] = turn_angle(ctx, 1);

    /* A B generator */
    /* ctx->gen_len = 2; */
    /* ctx->rot[0] = turn_angle(ctx, 1); */
    /* ctx->rot[1] = turn_angle(ctx, 1); */


    /* A^-4 B generator */
    /* ctx->gen_len = 2; */
    /* ctx->rot[0] = turn_angle(ctx, -4); */
    /* ctx->rot[1] = turn_angle(ctx, 1); */

    /* n=7 generator that seems to jumble */
    /* ctx->gen_len = 4; */

    /* ctx->rot[0] = turn_angle(ctx, 1); */
    /* ctx->rot[1] = turn_angle(ctx, 2); */
    /* ctx->rot[2] = turn_angle(ctx, 3); */
    /* ctx->rot[3] = turn_angle(ctx, 4); */

    /* n=7 experiment */
    /* ctx->gen_len = 4; */

    /* ctx->rot[0] = turn_angle(ctx, 2); */
    /* ctx->rot[1] = turn_angle(ctx, 2); */
    /* ctx->rot[2] = turn_angle(ctx, 3); */
    /* ctx->rot[3] = turn_angle(ctx, 3); */

    fprintf(stderr, "== RENDER SETTINGS ==\n");
    fprintf(stderr, "N: %d\n", ctx->n);
    fprintf(stderr, "R: %.15f\n", (double)ctx->r);
    fprintf(stderr, "Wedge-Only: %d\n", ctx->wedge_only);
    fprintf(stderr, "Box (%.5f, %.5f) to (%.5f, %.5f)\n",
            ctx->xmin, ctx->ymin,
            ctx->xmax, ctx->ymax);
    fprintf(stderr, "Image size %dx%d\n", ctx->img_w, ctx->img_h);
    fprintf(stderr, "Float precision set to %d, sizeof(FLOAT_T) is %lu, mantissa digits is %d\n", FPREC, sizeof(FLOAT_T), MANTDIG_L);



    /*test_xy_point(ctx);*/

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

    /* aa image */
    /*image_aa_sobel(ctx);*/

    fprintf(stderr, "Highest order seen: %lu\n", ctx->highest_order);

    ctx_to_png(ctx);


    free(ctx->grid);
    free(ctx);

    return 0;
}

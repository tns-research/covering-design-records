/* Classic post-optimisation baseline for (v,k,t) covering designs, no LLM.
 *
 * Deletes the block covering the fewest t-subsets alone, then runs simulated
 * annealing (Nurmela-Ostergard style moves) on the number of uncovered
 * t-subsets: pick an uncovered T, pick a block B meeting T in t-1 points, swap
 * the point of T missing from B for a point of B outside T. Repeats while it
 * succeeds and time remains. Prints the best complete design found.
 *
 * Usage: cd_postop v k t seconds seed < blocks.txt > out.txt
 *   blocks.txt: one block per line, 1-based points separated by spaces/commas.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int V, K, T;
static uint64_t binom[130][12];
static int nb, cap;
static int *blk;           /* nb x K, sorted ascending, 0-based */
static uint16_t *cnt;      /* per t-subset coverage count */
static long nunc;
static uint32_t *upos;     /* position in ulist, per t-subset */
static uint64_t *ulist;    /* uncovered t-subsets */

static uint64_t rng = 88172645463325252ULL;
static inline uint64_t rnd(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

static inline uint64_t rank(const int *s) { /* sorted ascending */
    uint64_t r = 0;
    for (int i = 0; i < T; i++) r += binom[s[i]][i + 1];
    return r;
}

/* apply f to every t-subset of block b (k points) */
static void each_sub(const int *b, int add) {
    int idx[12], s[12];
    for (int i = 0; i < T; i++) idx[i] = i;
    for (;;) {
        for (int i = 0; i < T; i++) s[i] = b[idx[i]];
        uint64_t r = rank(s);
        if (add > 0) {
            if (cnt[r]++ == 0) { uint64_t last = ulist[--nunc]; ulist[upos[r]] = last; upos[last] = upos[r]; }
        } else {
            if (--cnt[r] == 0) { ulist[nunc] = r; upos[r] = (uint32_t)nunc; nunc++; }
        }
        int i = T - 1;
        while (i >= 0 && idx[i] == K - T + i) i--;
        if (i < 0) break;
        idx[i]++;
        for (int j = i + 1; j < T; j++) idx[j] = idx[j - 1] + 1;
    }
}

static long uniq_of(const int *b) {
    int idx[12], s[12]; long u = 0;
    for (int i = 0; i < T; i++) idx[i] = i;
    for (;;) {
        for (int i = 0; i < T; i++) s[i] = b[idx[i]];
        if (cnt[rank(s)] == 1) u++;
        int i = T - 1;
        while (i >= 0 && idx[i] == K - T + i) i--;
        if (i < 0) break;
        idx[i]++;
        for (int j = i + 1; j < T; j++) idx[j] = idx[j - 1] + 1;
    }
    return u;
}

static void unrank(uint64_t r, int *s) {
    for (int i = T; i >= 1; i--) {
        int x = i - 1;
        while (binom[x + 1][i] <= r) x++;
        s[i - 1] = x;
        r -= binom[x][i];
    }
}

static void sort_block(int *b) {
    for (int i = 1; i < K; i++) { int x = b[i], j = i - 1; while (j >= 0 && b[j] > x) { b[j + 1] = b[j]; j--; } b[j + 1] = x; }
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: cd_postop v k t seconds seed < blocks\n"); return 2; }
    V = atoi(argv[1]); K = atoi(argv[2]); T = atoi(argv[3]);
    double secs = atof(argv[4]); rng ^= (uint64_t)atoll(argv[5]) * 0x9E3779B97F4A7C15ULL;
    for (int n = 0; n < 130; n++) { binom[n][0] = 1; for (int j = 1; j < 12; j++) binom[n][j] = n ? binom[n - 1][j - 1] + binom[n - 1][j] : 0; }
    uint64_t NT = binom[V][T];
    if (NT > 60000000ULL) { fprintf(stderr, "too many t-subsets\n"); return 2; }
    cnt = calloc(NT, sizeof(uint16_t)); nunc = (long)NT;
    upos = malloc(sizeof(uint32_t) * NT); ulist = malloc(sizeof(uint64_t) * NT);
    for (uint64_t i = 0; i < NT; i++) { ulist[i] = i; upos[i] = (uint32_t)i; }
    cap = 1024; blk = malloc(sizeof(int) * K * cap); nb = 0;
    int x, pos = 0;
    while (scanf("%d", &x) == 1) {
        if (nb == cap) { cap *= 2; blk = realloc(blk, sizeof(int) * K * cap); }
        blk[nb * K + pos++] = x - 1;
        if (pos == K) { sort_block(&blk[nb * K]); pos = 0; nb++; }
    }
    for (int i = 0; i < nb; i++) each_sub(&blk[i * K], 1);
    if (nunc) { fprintf(stderr, "input is not a covering (%ld uncovered)\n", nunc); return 1; }
    fprintf(stderr, "read %d blocks, C(%d,%d,%d)\n", nb, V, K, T);
    int best_nb = nb; int *best = malloc(sizeof(int) * K * nb); memcpy(best, blk, sizeof(int) * K * nb);
    clock_t t0 = clock();
    int s[12], nbk[25];
    while ((double)(clock() - t0) / CLOCKS_PER_SEC < secs) {
        /* delete weakest block */
        long bu = -1; int bi = 0;
        for (int i = 0; i < nb; i++) { long u = uniq_of(&blk[i * K]); if (bu < 0 || u < bu) { bu = u; bi = i; } }
        each_sub(&blk[bi * K], -1);
        memcpy(&blk[bi * K], &blk[(nb - 1) * K], sizeof(int) * K); nb--;
        double temp = 1.0; long it = 0;
        while (nunc > 0 && (double)(clock() - t0) / CLOCKS_PER_SEC < secs) {
            /* random uncovered t-subset */
            uint64_t r = ulist[rnd() % nunc];
            unrank(r, s);
            /* random block meeting s in t-1 points */
            int tries = 0, b = -1, out = -1, in = -1;
            while (tries++ < 4 * nb) {
                int c = rnd() % nb; int *B = &blk[c * K]; int m = 0, miss = -1;
                for (int i = 0, j = 0; i < T; i++) { while (j < K && B[j] < s[i]) j++; if (j < K && B[j] == s[i]) m++; else miss = s[i]; }
                if (m != T - 1) continue;
                int cand[25], nc = 0;
                for (int j = 0; j < K; j++) { int inS = 0; for (int i = 0; i < T; i++) if (s[i] == B[j]) inS = 1; if (!inS) cand[nc++] = B[j]; }
                b = c; in = miss; out = cand[rnd() % nc]; break;
            }
            if (b < 0) continue;
            int *B = &blk[b * K];
            long before = nunc;
            each_sub(B, -1);
            memcpy(nbk, B, sizeof(int) * K);
            for (int j = 0; j < K; j++) if (nbk[j] == out) nbk[j] = in;
            sort_block(nbk);
            each_sub(nbk, 1);
            long d = nunc - before;
            if (d <= 0 || urand() < exp(-d / temp)) memcpy(B, nbk, sizeof(int) * K);
            else { each_sub(nbk, -1); each_sub(B, 1); }
            if (++it % 1000 == 0) temp *= 0.999;
            if (temp < 0.05) temp = 1.0;
        }
        if (nunc == 0) {
            best_nb = nb; best = realloc(best, sizeof(int) * K * nb); memcpy(best, blk, sizeof(int) * K * nb);
            fprintf(stderr, "found %d blocks at %.0fs\n", nb, (double)(clock() - t0) / CLOCKS_PER_SEC);
        } else break;
    }
    for (int i = 0; i < best_nb; i++) { for (int j = 0; j < K; j++) printf(j ? " %d" : "%d", best[i * K + j] + 1); printf("\n"); }
    fprintf(stderr, "best %d\n", best_nb);
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAX_LEN 128

/* levenshtein with caller-supplied row buffers (no malloc per call) */
static int levenshtein(const unsigned char *a, size_t la,
                       const unsigned char *b, size_t lb,
                       size_t *prev, size_t *curr)
{
    for (size_t j = 0; j <= lb; j++) prev[j] = j;

    for (size_t i = 1; i <= la; i++) {
        curr[0] = i;
        for (size_t j = 1; j <= lb; j++) {
            if (a[i-1] == b[j-1]) {
                curr[j] = prev[j-1];
            } else {
                size_t sub = prev[j-1] + 1;
                size_t del = prev[j]   + 1;
                size_t ins = curr[j-1] + 1;
                curr[j] = sub < del ? (sub < ins ? sub : ins)
                                    : (del < ins ? del : ins);
            }
        }
        /* swap rows */
        size_t *tmp = prev; prev = curr; curr = tmp;
    }
    return (int)prev[lb];
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "usage: %s file\n", argv[0]);
        return 1;
    }

    FILE *fp = fopen(argv[1], "r");
    if (!fp) {
        perror("fopen");
        return 1;
    }

    /* --- load sequences --- */
    size_t capacity = 1024;
    size_t n        = 0;

    unsigned char *seqs = aligned_alloc(64, capacity * MAX_LEN);
    size_t        *lens = malloc(capacity * sizeof(size_t));
    if (!seqs || !lens) {
        perror("alloc");
        return 1;
    }

    char buf[MAX_LEN];
    while (fgets(buf, sizeof(buf), fp)) {
        buf[strcspn(buf, "\n")] = 0;
        size_t l = strlen(buf);
        if (l == 0) continue;
        if (l >= MAX_LEN) {
            fprintf(stderr,
                    "sequence %zu exceeds MAX_LEN=%d, skipping\n",
                    n, MAX_LEN);
            continue;
        }

        if (n >= capacity) {
            capacity *= 2;
            seqs = realloc(seqs, capacity * MAX_LEN);
            lens = realloc(lens, capacity * sizeof(size_t));
            if (!seqs || !lens) {
                perror("realloc");
                return 1;
            }
        }

        memcpy(seqs + n * MAX_LEN, buf, l);
        lens[n] = l;
        n++;
    }
    fclose(fp);

    if (n == 0) {
        fprintf(stderr, "no sequences loaded\n");
        free(seqs); free(lens);
        return 1;
    }

    fprintf(stderr, "loaded %zu sequences\n", n);

    /* --- allocate thread-local row buffers --- */
    int n_threads = omp_get_max_threads();
    /* 2 rows per thread, each MAX_LEN+1 elements */
    size_t *row_buffers = malloc(
        2 * (size_t)n_threads * (MAX_LEN + 1) * sizeof(size_t)
    );
    if (!row_buffers) {
        perror("alloc row_buffers");
        return 1;
    }

    /* --- parallel pairwise Levenshtein --- */
    double   global_sum    = 0.0;
    double   global_sumsq  = 0.0;
    uint64_t global_pairs  = 0;

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        size_t *t_prev = row_buffers + (2 * tid)     * (MAX_LEN + 1);
        size_t *t_curr = row_buffers + (2 * tid + 1) * (MAX_LEN + 1);

        double   local_sum   = 0.0;
        double   local_sumsq = 0.0;
        uint64_t local_pairs = 0;

        #pragma omp for schedule(dynamic, 32)
        for (size_t i = 0; i < n; i++) {
            unsigned char *a = seqs + i * MAX_LEN;
            size_t la = lens[i];

            for (size_t j = i + 1; j < n; j++) {
                unsigned char *b = seqs + j * MAX_LEN;
                size_t lb = lens[j];
                size_t max_len = (la > lb) ? la : lb;

                int dist = levenshtein(a, la, b, lb, t_prev, t_curr);
                double p = (double)dist / (double)max_len;
                local_sum   += p;
                local_sumsq += p * p;
                local_pairs++;
            }
        }

        #pragma omp atomic
        global_sum += local_sum;
        #pragma omp atomic
        global_sumsq += local_sumsq;
        #pragma omp atomic
        global_pairs += local_pairs;
    }

    /* --- compute and print statistics --- */
    double mean     = global_sum / (double)global_pairs;
    double variance = (global_sumsq / (double)global_pairs) - (mean * mean);

    fprintf(stderr, "pairs     : %lu\n",    global_pairs);
    fprintf(stderr, "mean      : %.12f\n",  mean);
    fprintf(stderr, "variance  : %.12f\n",  variance);
    printf("%zu,%lu,%.12f,%.12f\n",
           n, global_pairs, mean, variance);

    free(seqs);
    free(lens);
    free(row_buffers);
    return 0;
}

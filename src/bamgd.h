#pragma once
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>
#include <inttypes.h>
#include <sys/ioctl.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include "ketopt.h"
#include "khashl.h"
#include "cgranges.h"
#include "thpool.h"
#include "kde.h"
#include "version.h"

extern const char *__progname;

#define CHUNK 0xFFFF
#define MAXDP 0xFFFFF
#define MINQL 35
#define MINMQ 10
#define WIN 1000

#define ARR "\e[90m\xE2\x97\x82\e[0m"
#define BUL "\e[90m\xE2\x97\x8F\e[0m"
#define INF "\E[1;34m\xE2\x84\xb9\E[0;0m"
#define ERR "\e[31m\xE2\x9C\x97\e[0m"
#define PP fprintf(stderr, "%s\t%d\t<%s>\n", __FILE__, __LINE__, __func__);

// set of strings
KHASHL_CSET_INIT(KH_LOCAL, ss_t, ss, const char *, kh_hash_str, kh_eq_str)

#define error(msg, ...) do {                                 \
    fputs("\e[31m\xE2\x9C\x97\e[0m ", stderr);               \
    fprintf(stderr, msg, ##__VA_ARGS__);                     \
    if (errno) fprintf(stderr, ": %s", strerror(errno));     \
    fflush(stderr);                                          \
    exit(EXIT_FAILURE);                                      \
} while (0)

#define PPT do \
{ \
    char buf[80]; \
    time_t now = time(0); \
    strftime(buf, sizeof(buf), "\e[34m%D %X\e[0m", localtime(&now)); \
    fprintf(stderr, "%s %s %s\t%d\t<%s>\n", INF, buf, __FILE__, __LINE__, __func__); \
} while (0);

// argument struct
typedef struct
{
	char *in, *out, *sub, *ref, *ctg;
	bool log, dup;
} arg_t;

static ko_longopt_t long_options[] = {
	{ "in",        ko_required_argument, 'i' },
	{ "out",       ko_required_argument, 'o' },
	{ "sub",       ko_required_argument, 's' },
	{ "ref",       ko_required_argument, 'r' },
	{ "ctg",       ko_required_argument, 'c' },
	{ "log",       ko_no_argument, 'l' },
	{ "dup",       ko_no_argument, 'd' },
	{ "help",      ko_no_argument, 'h' },
	{ "version",   ko_no_argument, 'v' },
	{ NULL, 0, 0 }
};

typedef struct
{
	char *in;
	char *ctg;
	bool dedup;
	cgranges_t *cr;
	double **dep;
} op_t;

typedef struct // auxiliary data structure
{
	samFile *fp;     // the file handle
	bam_hdr_t *hdr;  // the file header
	hts_itr_t *iter; // NULL if a region not specified
	int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

void ss_ins(ss_t *h, const char *a);
void ss_des(ss_t *h);
int read_bam(void *data, bam1_t *b);
int read_ddp_bam(void *data, bam1_t *b);
bool bam_has_dup(const char *fn);
void bam_get_ref(bam_hdr_t *h, char *ref);
void chk_fai(const arg_t *arg);
void prs_arg(int argc, char **argv, arg_t *arg);
void ld_cr(const faidx_t *fai, const char *ctg, cgranges_t *cr);
double seq_gc(const char *seq);
void cr_gc(const cgranges_t *cr, const faidx_t *fai, double *gc);
void ld_dp(void *op);
void draw_axis(cairo_t *cr, uint32_t md, const char *ctg, uint32_t n_targets,
		uint64_t gl);
int is_gzip(const char *fn);
bool ends_with(const char *str, const char *sfx);
int strlen_wo_esc(const char *str);
void horiz(int n, bool bold);
void usage();

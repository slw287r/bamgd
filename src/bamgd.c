#include "bamgd.h"

int main(int argc, char *argv[])
{
	if (argc == 1)
		usage();
	setenv("FONTCONFIG_PATH", "/etc/fonts", 1);
	arg_t *arg = calloc(1, sizeof(arg_t));
	prs_arg(argc, argv, arg);
	samFile *fp = sam_open(arg->in, "r");
	if (!fp)
		error("Error: failed to read input bam [%s]\n", arg->in);
	int i, j = 0, k = 0;
	cgranges_t *cr = cr_init();
	faidx_t *fai = fai_load(arg->ref);
	ld_cr(fai, arg->ctg, cr);
	dm_t dm = {100, 0, MAXDP, 0};
	double *gc = calloc(cr->n_r, sizeof(double));
	double *dp = calloc(cr->n_r, sizeof(double));
	/* dbg cr
	for (i = 0; i < cr->n_r; ++i)
	{
		cr_intv_t *r = &cr->r[i];
		// r->x and r->y are invalid after cr_index, use cr_st and cr_en instead
		printf("%s\t%d\t%d\t%d\t%d\t%d\n", cr->ctg[r->x>>32].name, cr_st(r), cr_en(r), (int32_t)r->x, r->y, r->label);
	}
	*/
	cr_gc(cr, fai, gc);
	/* debug gc
	for (i = 0; i < cr->n_r; ++i)
		printf("%f\n", gc[i]);
	*/
	cr_index(cr); // index cr for overlapping, no merge is required for ref fai
	ld_dp(arg, cr, dp);
	// shrink dp and gc
	for (i = 0; i < cr->n_r; ++i)
	{
		if (dp[i])
		{
			gc[j] = gc[i];
			dp[j++] = dp[i];
		}
	}
	gc = realloc(gc, j * sizeof(double));
	dp = realloc(dp, j * sizeof(double));
	double *dp_cp = calloc(j, sizeof(double));
	memcpy(dp_cp, dp, j * sizeof(double));
	gsl_sort(dp_cp, 1, j);
	double lower = gsl_stats_quantile_from_sorted_data (dp_cp, 1, j, 0.05);
	double upper = gsl_stats_quantile_from_sorted_data (dp_cp, 1, j, 0.95);
	free(dp_cp);
	//quantile filter of dp
	for (i = 0; i < j; ++i)
	{
		if (dp[i] >= lower && dp[i] <= upper)
		{
			dp[k] = dp[i];
			gc[k++] = gc[i];
		}
	}
	dp = realloc(dp, k * sizeof(double));
	gc = realloc(gc, k * sizeof(double));
	if (arg->log)
		for (i = 0; i < k; ++i)
			dp[i] = log10(dp[i]);
	for (i = 0; i < k; ++i)
	{
		dm.xmin = fmin(dm.xmin, gc[i]);
		dm.xmax = fmax(dm.xmax, gc[i]);
		dm.ymin = fmin(dm.ymin, dp[i]);
		dm.ymax = fmax(dm.ymax, dp[i]);
	}
	dm.ymin = arg->log ? 0.0 : 1.0;
	//printf("%f\t%f\t%f\t%f\n", dm.xmin, dm.xmax, dm.ymin, dm.ymax);
	kde_plot(gc, dp, k, &dm, arg->log, arg->out);
	/* dbg gc~dp
	for (i = 0; i < k; ++i)
		printf("%f\t%f\n", gc[i], dp[i]);
	*/
	cr_destroy(cr);
	free(dp);
	free(gc);
	free(arg);
	return 0;
}

int read_bam(void *data, bam1_t *b)
{
	aux_t *aux = (aux_t*)data;
	int ret;
	while (true)
	{
		ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
		if (ret < 0)
			break;
		if (b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
			continue;
		/*
		if ((int)b->core.qual < aux->min_mapQ)
			continue;
		if (aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len)
			continue;
		*/
		break;
	}
	return ret;
}

void ld_dp(const arg_t *arg, const cgranges_t *cr, double *dp)
{
	int tid, pos, beg = 0, end = INT_MAX, n_plp, i;
	int64_t m = 0, *b = 0, n_b = 0;
	uint64_t *cov = calloc(cr->n_r, sizeof(uint64_t));
	aux_t *data = calloc(1, sizeof(aux_t));
	data->fp = hts_open(arg->in, "r");
	data->hdr = sam_hdr_read(data->fp);
	hts_idx_t *idx = sam_index_load(data->fp, arg->in);
	if (arg->ctg)
	{
		data->iter = sam_itr_querys(idx, data->hdr, arg->ctg);
		beg = data->iter->beg;
		end = data->iter->end;
	}
	hts_idx_destroy(idx);
	bam_hdr_t *h = data->hdr;
	bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void**)&data);
	bam_mplp_set_maxcnt(mplp, MAXDP);
	const bam_pileup1_t *plp = calloc(1, sizeof(bam_pileup1_t));
	while (bam_mplp_auto(mplp, &tid, &pos, &n_plp, &plp) > 0)
	{
		if (pos < beg)
			continue;
		if (pos >= end || tid >= h->n_targets)
			break;
		if (!(n_b = cr_overlap(cr, sam_hdr_tid2name(h, tid), pos, pos, &b, &m)))
			continue;
		else
		{
			i = cr_label(cr, b[0]);
			dp[i] += n_plp;
			++cov[i];
		}
	}
	if (b)
		free(b);
	for (i = 0; i < cr->n_r; ++i)
		if (cov[i])
			dp[i] /= cov[i];
	free(cov);
	bam_hdr_destroy(data->hdr);
	if (data->fp)
		hts_close(data->fp);
	hts_itr_destroy(data->iter);
	free(data);
}

int is_gzip(const char *fn)
{
	char buf[2];
	FILE *fp;
	int gzip = 0;
	if ((fp = fopen(fn, "rb")) == NULL)
		error("[ERROR] Unable to open file: %s\n", fn);
	if (fread(buf, 1, 2, fp) == 2)
		if (((int)buf[0] == 0x1f) && ((int)(buf[1]&0xFF) == 0x8b))
			gzip = 1;
	fclose(fp);
	return gzip;
}

/*
 * @PG  ID:bwa       PN:bwa       VN:0.7.18  CL:bwa mem -at8 /path/to/ref.fa /path/to/fq.gz -R
 * @PG  ID:bwa-mem2  PN:bwa-mem2  VN:2.2.1a  CL:bwa-mem2 mem -v1 -zt8 -h64 /path/to/ref.fa /path/to/fq.gz
 */
void bam_get_ref(bam_hdr_t *h, char *ref)
{
	int i;
	kstring_t ks = {0, 0, 0};
	if (sam_hdr_find_line_id(h, "PG", "ID", "bwa-mem2", &ks) &&
			sam_hdr_find_line_id(h, "PG", "ID", "bwa", &ks))
		return;
	const char sfx[][16] = {".fasta.gz ", ".fasta ", ".fa.gz ", ".fa ", ".fna.gz ", ".fna "};
	char *p = NULL, *q = NULL;
	for (i = 0; i < 6; ++i)
	{
		if ((p = strstr(ks.s, sfx[i])))
		{
			q = p;
			while (!isspace(*q--));
			strncpy(ref, q + 2, p - q + strlen(sfx[i]) - 3);
			break;
		}
	}
	if (ks.m) free(ks.s);
}

void ss_ins(ss_t *h, const char *a)
{
	khint_t k;
	int absent;
	if ((k = ss_get(h, a)) == kh_end(h))
		k = ss_put(h, strdup(a), &absent);
}

void ss_des(ss_t *h)
{
	khint_t k;
	kh_foreach(h, k)
		free((char *)kh_key(h, k));
	ss_destroy(h);
}

bool ends_with(const char *str, const char *sfx)
{
	int ret = 0;
	int str_len = strlen(str);
	int sfx_len = strlen(sfx);
	if ((str_len >= sfx_len) && (0 == strcasecmp(str + (str_len-sfx_len), sfx)))
		ret = 1;
	return ret;
}

void chk_fai(const arg_t *arg)
{
	if (access(arg->ref, R_OK))
		error("Error: specified reference [%s] is inaccessible!\n", arg->ref);
	faidx_t *fai = NULL;
	if (!(fai = fai_load(arg->ref)))
		error("Error loading reference index of [%s]\n", arg->ref);
	if (arg->ctg && !faidx_has_seq(fai, arg->ctg))
	{
		fai_destroy(fai);
		error("Specified contig [%s] not found in reference [%s]!\n",
				arg->ctg, arg->ref);
	}
	samFile *fp = sam_open(arg->in, "r");
	bam_hdr_t *hdr = sam_hdr_read(fp);
	if (!hdr)
		error("Error reading input bam header from [%s]!\n", arg->in);
	if (arg->ctg && sam_hdr_name2tid(hdr, arg->ctg) < 0)
		error("Specified contig [%s] is not found in bam!\n", arg->ctg);
	khint_t k;
	int i, n = 0;
	ss_t *s = ss_init();
	for (i = 0; i < faidx_nseq(fai); ++i)
	{
		const char *ctg = faidx_iseq(fai, i);
		if (sam_hdr_name2tid(hdr, ctg) < 0)
			ss_ins(s, ctg);
	}
	fai_destroy(fai);
	if ((n = kh_size(s)) >= 100)
	{
		ss_des(s);
		error("Too much reference sequences (%d in total) not found in bam!", n);
	}
	else if (n > 0)
	{
		fprintf(stderr, "The following %d reference sequence%s not found in bam:\n",
				n, n > 1 ? "s are" : " is");
		kh_foreach(s, k)
		{
			fputs(kh_key(s, k), stderr);
			fputc(' ', stderr);
		}
		fputc('\n', stderr);
		ss_des(s);
		error("Reference error!\n");
	}
	ss_des(s);
	bam_hdr_destroy(hdr);
	hts_close(fp);
}

void prs_arg(int argc, char **argv, arg_t *arg)
{
	int c = 0;
	ketopt_t opt = KETOPT_INIT;
	const char *opt_str = "i:o:s:r:c:lhv";
	while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0)
	{
		switch (c)
		{
			case 'i': arg->in = opt.arg; break;
			case 'o': arg->out = opt.arg; break;
			case 's': arg->sub = opt.arg; break;
			case 'r': arg->ref = opt.arg; break;
			case 'c': arg->ctg = opt.arg; break;
			case 'l': arg->log = true; break;
			case 'h': usage(); break;
			case 'v':
				if (strlen(BRANCH_COMMIT))
					printf("%s (%s)\n", VERSION, BRANCH_COMMIT);
				else
					puts(VERSION);
				exit(EXIT_SUCCESS);
			case '?':
				printf("Invalid option: [%c]\n", opt.opt); exit(EXIT_SUCCESS);
			default:
				printf("Invalid option: [%s]", opt.arg); exit(EXIT_SUCCESS);
		}
	}
	if (!arg->in || access(arg->in, R_OK))
		error("Error: input bam is unspecified or inaccessible!\n");
	if (!(ends_with(arg->in, ".bam") && is_gzip(arg->in)))
		error("Oops! only bam input is supported. Invalid bam file [%s]\n", arg->in);
	char bai[PATH_MAX];
	snprintf(bai, PATH_MAX, "%s.bai", arg->in);
	if (access(bai, R_OK))
		error("Error: bam's index file (.bai) is required, please use samtools sort and index to create it.\n");
	if (!arg->ref)
	{
		static char ref[PATH_MAX] = {'\0'};
		samFile *fp = sam_open(arg->in, "r");
		bam_hdr_t *hdr = sam_hdr_read(fp);
		if (!hdr)
			error("Error reading input bam header!\n");
		if (arg->ctg && sam_hdr_name2tid(hdr, arg->ctg) < 0)
			error("Specified contig [%s] is not found in bam!\n", arg->ctg);
		bam_get_ref(hdr, ref);
		bam_hdr_destroy(hdr);
		hts_close(fp);
		if (!strlen(ref))
			error("Failed getting unspecified reference from bam file!\n");
		if (access(ref, R_OK))
			error("Reference [%s] in bam is inaccessible!\n", ref);
		arg->ref = ref;
	}
	chk_fai(arg);
	if (!arg->out)
	{
		static char png[PATH_MAX];
		char *p = strrchr(arg->in, '/');
		if (p)
			strncpy(png, p + 1, PATH_MAX);
		else
			strncpy(png, arg->in, PATH_MAX);
		p = strrchr(png, '.');
		strncpy(p + 1, "png", 3);
		arg->out = png;
	}
	else if (!ends_with(arg->out, ".png"))
		error("Only png output is supported, got [%s]\n", arg->out);
}

void ld_cr(const faidx_t *fai, const char *ctg, cgranges_t *cr)
{
	int idx = 0;
	uint64_t i, j, gl = 0, cl = 0;
	if (ctg)
		gl = faidx_seq_len64(fai, ctg);
	else
		for (i = 0; i < faidx_nseq(fai); ++i)
			gl += faidx_seq_len64(fai, faidx_iseq(fai, i));
	if (!gl)
		error("Error getting genome size!\n");
	int win = gl > WIN ? WIN : pow(10, (int)log10(gl) - 1);
	if (ctg)
		for (i = 0; i < gl; i += win)
			cr_add(cr, ctg, i, (uint64_t)fmin(i + win, gl), idx++);
	else
	{
		for (i = 0; i < faidx_nseq(fai); ++i)
		{
			const char *c1 = faidx_iseq(fai, i);
			cl = faidx_seq_len64(fai, c1);
			for (j = 0; j < cl; j += win)
				cr_add(cr, c1, j, (uint64_t)fmin(j + win, cl), idx++);
		}
	}
}

double seq_gc(const char *seq)
{
	uint64_t bases[5] = {0ull};
	const char *p = seq;
	while (*p++)
		++bases[seq_nt16_int[seq_nt16_table[(int)*p]]];
	return 100.0 * (bases[1] + bases[2]) / (bases[0] + bases[1] + bases[2] + bases[3]);
}

void cr_gc(const cgranges_t *cr, const faidx_t *fai, double *gc)
{
	int i, l;
	char *seq;
	for (i = 0; i < cr->n_r; ++i)
	{
		cr_intv_t *r = &cr->r[i];
		seq = faidx_fetch_seq(fai, cr->ctg[r->x>>32].name, (int32_t)r->x, r->y, &l);
		gc[i] = l ? seq_gc(seq) : 0;
		free(seq);
	}
}

int strlen_wo_esc(const char *str)
{
	int length = 0;
	while (*str != '\0')
	{
		// Check if this is the start of an ANSI escape sequence (ESC [ ... m)
		if (*str == '\033' && *(str + 1) == '[')
		{
			// Move past \033[
			str += 2;
			// Skip until we reach 'm', which ends the ANSI sequence
			while (*str != '\0' && *str != 'm')
				str++;
			// Move past 'm' if we found it
			if (*str == 'm')
				str++;
			continue;
		}
		// Count visible characters
		length++;
		str++;
	}
	return length;
}

void horiz(int n, bool bold)
{
	struct winsize w;
	ioctl(0, TIOCGWINSZ, &w);
	n = (w.ws_col >= n) ? n : w.ws_col;
	printf("\e[%s90m", bold ? "1;" : "");
	while (n--) fputs("\xe2\x94\x80", stdout);
	printf("\e[%s0m\n", bold ? "0;" : "");
}

void usage()
{
	int w = 58;
	horiz(w, true);
	char title[] = "\e[1mPlot genome GC vs sequencing depth in bam\e[0m";
	int title_len = strlen_wo_esc(title);
	printf("%*.*s\n", (int)((w - title_len) / 2 + strlen(title)), (int)strlen(title), title);
	horiz(w, false);
	printf("%s \e[1mUsage\e[0m: \e[1;31m%s\e[0;0m \e[1;90m[options]\e[0;0m --in <bam> --out <png>\n", BUL, __progname);
	putchar('\n');
	puts(BUL " \e[1mOptions\e[0m:");
	puts("  -i, --in  \e[3mFILE\e[0m   Input BAM file with bai index");
	puts("  -o, --out \e[3mSTR\e[0m    Output GC~depth plot png \e[90m[${prefix}.png]\e[0m");
	puts("  -r, --ref \e[3mSTR\e[0m    Reference fasta used in bam \e[90m[auto]\e[0m");
	puts("  -c, --ctg \e[3mSTR\e[0m    Restrict plot to this contig \e[90m[none]\e[0m");
	puts("  -l, --log \e[3mBOOL\e[0m   Depth in logscale \e[90m[false]\e[0m");
	putchar('\n');
	puts("  -h               Show help message");
	puts("  -v               Display program version");
	putchar('\n');
	puts(BUL " \e[1mContact\e[0m: \e[4mmeta@geneplus.cn\e[0m for support and bug report");
	horiz(w, true);
	exit(1);
}

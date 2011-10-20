#include "output_vcf.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/* vcf header fields */
static const char *_vcf_format_string = "GT:DP:DP4:GQ:VAQ:BQ:MQ:SS:SSC";
const static struct {
    const char *id;
    const char *number;
    const char *type;
    const char *description;
} _vcf_format_fields[] = {
    { "GT", "1", "String", "Genotype" },
    { "DP", "1", "Integer", "Total read depth" },
    { "DP4", "4", "Integer", "# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases" },
    { "GQ", "1", "Integer", "Genotype quality" },
    { "VAQ", "1", "Integer", "Variant allele quality" },
    { "BQ", "G", "Integer", "Average base quality" },
    { "MQ", "G", "Integer", "Average mapping quality" },
    { "SS", "1", "Integer",
        "Variant status relative to non-adjacent Normal, 0=wildtype,1=germline,2=somatic,3=LOH,4=unknown" },
    { "SSC", "1", "Integer", "Somatic Score" },
};
const static uint32_t _n_vcf_format_fields = sizeof(_vcf_format_fields)/sizeof(_vcf_format_fields[0]);

static void output_vcf_intN(FILE *fh, const uint32_t *values, uint32_t n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (i > 0)
            fputc(',', fh);
        fprintf(fh, "%d", values[i]);
    }
}

static void output_vcf_int4_masked(FILE *fh, const uint32_t values[4], uint32_t mask) {
    uint32_t i;
    uint32_t c = 0;
    for (i = 0; i < 4; ++i) {
        uint32_t val = 1 << i;
        if (mask & val) {
            if (++c > 1)
                fputc(',', fh);
            fprintf(fh, "%d", values[i]);
        }
    }
}

static void output_vcf_gt(FILE *fh, uint32_t ref_base, uint32_t alts, uint32_t gt) {
    uint32_t allele_idx = 0;
    uint32_t out_count = 0;
    uint32_t allele_count = 0;
    uint32_t i;

    for (i = 0; i < 4; ++i) {
        int value = 1 << i;
        if (gt & value)
            ++allele_count;
    }

    if (gt & ref_base) {
        if (allele_count == 1) {
            fprintf(fh, "0/0");
            return;
        }

        fprintf(fh, "0");
        ++out_count;
    }
    gt &= ~ref_base;

    for (i = 0; i < 4; ++i) {
        int value = 1 << i;
        if (alts & value)
            ++allele_idx;

        if (gt & value) {
            if (allele_count == 1) {
                fprintf(fh, "%d/%d", allele_idx, allele_idx);
                return;
            }
            if (out_count > 0)
                fputc('/', fh);
            fprintf(fh, "%d", allele_idx);
            ++out_count;
        }
    }
}

void output_vcf_sample(FILE *fh, int ref_base4, int alts, const sample_data_t *s) {
    /* GT */
    output_vcf_gt(fh, ref_base4, alts, s->genotype);

    /* DP */
    fprintf(fh, ":%d:", s->dqstats.total_depth);

    /* DP4 */
    output_vcf_intN(fh, s->dqstats.dp4, 4);

    /* GQ, VAQ */
    fprintf(fh, ":%d:%d:", s->consensus_quality, s->variant_allele_quality);

    /* BQ, MQ */
    output_vcf_int4_masked(fh, s->dqstats.mean_baseQ, s->genotype);
    fputc(':', fh);
    output_vcf_int4_masked(fh, s->dqstats.mean_mapQ, s->genotype);
    fputc(':', fh);

    /* SS */
    fputc(s->is_somatic ? '2' : '.', fh);
    fputc(':', fh);

    /* SSC */
    if (s->somatic_score >= 0) {
        fprintf(fh, "%d", s->somatic_score);
    } else {
        fputc('.', fh);
    }
}

void output_vcf(FILE *fh, const sniper_output_t *p) {
    uint32_t n_alt = 0;
    uint32_t alts = (p->tumor.genotype | p->normal.genotype) & ~p->ref_base4;
    uint32_t i;

    /* CHROM, POS, ID, REF */
    fprintf(fh, "%s\t%d\t.\t%c\t", 
        p->seq_name,
        p->pos + 1,
        p->ref_base
    );

    /* ALT */
    for (i = 0; i < 4; ++i) {
        int value = 1 << i;
        if (alts & value) {
            if (++n_alt > 1)
                fputc(',', fh);
            fputc(bam_nt16_rev_table[value], fh);
        }
    }
    if (n_alt == 0)
        fputc('.', fh);
    fputc('\t', fh);

    /* QUAL and FILTER fields */
    fprintf(fh, ".\t.\t");

    /* INFO, FORMAT */
    fprintf(fh, ".\t%s\t", _vcf_format_string);

    output_vcf_sample(fh, p->ref_base4, alts, &p->normal);
    fputc('\t', fh);
    output_vcf_sample(fh, p->ref_base4, alts, &p->tumor);
    fputc('\n', fh);
}

void output_vcf_header(FILE *fh, const header_data_t *h) {
    char filedate[255];
    time_t now = time(NULL);
    uint32_t i;
    strftime(filedate, sizeof(filedate), "%Y%m%d", localtime(&now));
    fprintf(fh, "##fileformat=VCFv4.1\n");
    fprintf(fh, "##fileDate=%s\n", filedate);
    fprintf(fh, "##reference=file://%s\n", h->refseq);
    for (i = 0; i < _n_vcf_format_fields; ++i) {
        fprintf(fh, "##FORMAT=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n",
            _vcf_format_fields[i].id,
            _vcf_format_fields[i].number,
            _vcf_format_fields[i].type,
            _vcf_format_fields[i].description
            );
    }
    fprintf(fh, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s	%s\n",
        h->normal_sample_id,
        h->tumor_sample_id);
}

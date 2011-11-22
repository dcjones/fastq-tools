
bin_PROGRAMS = fastq-grep fastq-kmers fastq-match fastq-uniq fastq-qual fastq-sample fastq-qualadj

fastq_common_src=common.h common.c
fastq_parse_src=parse.h parse.c
fastq_sw_src=sw.h sw.c
fastq_hash_src=hash.h hash.c
fastq_rng_src=rng.h rng.c

fastq_grep_SOURCES = fastq-grep.c $(fastq_common_src) $(fastq_parse_src)
fastq_grep_LDADD = $(PCRE_LIBS) -lz

fastq_kmers_SOURCES = fastq-kmers.c $(fastq_common_src) $(fastq_parse_src)
fastq_kmers_LDADD = -lz

fastq_match_SOURCES = fastq-match.c $(fastq_common_src) $(fastq_parse_src) $(fastq_sw_src)
fastq_match_LDADD   = -lz

fastq_uniq_SOURCES = fastq-uniq.c $(fastq_common_src) $(fastq_parse_src) $(fastq_hash_src)
fastq_uniq_LDADD   = -lz

fastq_qual_SOURCES = fastq-qual.c $(fastq_common_src) $(fastq_parse_src)
fastq_qual_LDADD   = -lz

fastq_sample_SOURCES = fastq-sample.c $(fastq_common_src) $(fastq_parse_src) $(fastq_rng_src)
fastq_sample_LDADD   = -lz

fastq_qualadj_SOURCES = fastq-qualadj.c $(fastq_common_src) $(fastq_parse_src)
fastq_qualadj_LDADD   = -lz

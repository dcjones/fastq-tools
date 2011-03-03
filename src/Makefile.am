
bin_PROGRAMS = fastq-grep fastq-kmers fastq-match

fastq_common_src=fastq-common.h fastq-common.c
fastq_parse_src=fastq-parse.h fastq-parse.c
fastq_sw_src=fastq-sw.h fastq-sw.c

fastq_grep_SOURCES = fastq-grep.c $(fastq_common_src) $(fastq_parse_src)
fastq_grep_LDADD = $(PCRE_LIBS) -lz

fastq_kmers_SOURCES = fastq-kmers.c $(fastq_common_src) $(fastq_parse_src)
fastq_kmers_LDADD = -lz

fastq_match_SOURCES = fastq-match.c $(fastq_common_src) $(fastq_parse_src) $(fastq_sw_src)
fastq_match_LDADD   = -lz



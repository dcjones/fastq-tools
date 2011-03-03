
SUBDIRS = swsse2

bin_PROGRAMS = fastq-grep fastq-kmers fastq-match

fastq_grep_SOURCES = kseq.h fastq-grep.c fastq-parse.c fastq-common.c
fastq_grep_LDADD = $(PCRE_LIBS) -lz

fastq_kmers_SOURCES = kseq.h fastq-kmers.c fastq-parse.c fastq-common.c
fastq_kmers_LDADD = -lz

fastq_match_SOURCES = kseq.h fastq-match.c fastq-parse.c fastq-common.c
fastq_match_LDADD   = swsse2/libswsse2.la -lz
fastq_match_CFLAGS  = -msse2


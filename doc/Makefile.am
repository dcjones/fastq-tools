
dist_man_MANS = fastq-grep.1 fastq-kmers.1 fastq-match.1 fastq-uniq.1 fastq-sample.1 fastq-sort.1

.1.html :
	groff -man -Thtml < $< \
		| perl -pe 's/<head>/<head>\n<link rel="stylesheet" href="stylesheets\/man.css">/' \
		> $@

html-local : fastq-grep.html fastq-kmers.html fastq-match.html fastq-uniq.html fastq-sample.html fastq-sort.html

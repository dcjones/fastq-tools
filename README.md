
fastq-tools
===========

This package provides a number of small and efficient programs to perform common
tasks with high throughput sequencing data in the FASTQ format. All of the
programs work with typical FASTQ files as well as gzipped FASTQ files.


index
-----

The following programs are provided. See the individual man pages for more
information.

* *fastq-sort* : sort fastq entries by various keys

* *fastq-grep* : match sequences against regular expressions

* *fastq-kmers* : count k-mer occurrences

* *fastq-match* : (smith-waterman) local sequence alignment

* *fastq-qual* : tabulate quality scores

* *fastq-sample* : randomly sample reads, with or without replacement

* *fastq-uniq* : count duplicate reads

* *fastq-qualadj* : adjust quality scores by a fixed offset


install
-------

On most systems, installation is as simple as `./configure && make install`.

If the source was obtained from the git repository, the included `./autogen.sh`
script must be run first to generate the `configure` script.

The only external dependencies are PCRE (http://www.pcre.org/) and zlib
(http://zlib.net/).

    $ sudo apt-get install libpcre3-dev
    $ sudo apt-get install zlibc
    $ git clone https://github.com/dcjones/fastq-tools.git
    $ cd fastq-tools
    $ ./autogen.sh
    $ ./configure
    $ make install
    $ cp -puv src/fastq-{grep,kmers,match,uniq,qual,sample,qualadj,sort,qscale} /usr/local/bin/
    $ cd ..
    $ rm fastq-tools -rf


contribute
----------

If you have written any small but useful programs to deal with FASTQ files,
please consider submitting them for inclusion in fastq-tools. Check out the
Github page (https://github.com/dcjones/fastq-tools) or send mail to the author
(dcjones@cs.washington.edu).


copying
-------

This package is provided under a permissive MIT-style license. In particular:

Copyright (C) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.



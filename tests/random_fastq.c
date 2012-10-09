
/*
    random_fastq
    ------------
    Generate random data in FASTQ format.
*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

static void print_help()
{
    printf(
"Usage: random_fastq [option]...\n"
"Generate an endless stream of random FASTQ data to standard out.\n\n"
"Options:\n"
"    TODO\n\n"
"Beware: the only purpose of this program is test quip.\n"
"No particular guarantees are made.\n\n");
}


/* Draw n samples from a categorial distribution of size k with cumulative
 * distribution given by cs and element us. The sample is stored in xs which
 * is assumed to be of appropriate length. */
void randcat(const char* us, const double* cs, size_t k, char* xs, size_t n)
{
    size_t i;
    double r;
    for (i = 0; i < n; ++i) {
        r = drand48();
        int a = 0, b = (int) k, c;
        do {
            c = (a + b) / 2;
            if (r <= cs[c]) b = c;
            else            a = c + 1;
        } while (a < b);

        xs[i] = us[a];
    }
}


int main(int argc, char* argv[])
{
    static struct option long_options[] =
    {
        {"min-length", required_argument, NULL, 'm'},
        {"max-length", required_argument, NULL, 'M'},
        {"length",     required_argument, NULL, 'l'},
        {"id-length",  required_argument, NULL, 'i'},
        {"help",       no_argument,       NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    size_t min_len = 100, max_len = 100;
    size_t id_len = 50;

    int opt, opt_idx;
    while (1) {
        opt = getopt_long(argc, argv, "h", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'm':
                min_len = (size_t) strtoul(optarg, NULL, 10);
                break;

            case 'M':
                max_len = (size_t) strtoul(optarg, NULL, 10);
                break;

            case 'l':
                min_len = max_len = (size_t) strtoul(optarg, NULL, 10);
                break;

            case 'i':
                id_len = (size_t) strtoul(optarg, NULL, 10);
                break;

            case 'h':
                print_help();
                return EXIT_SUCCESS;

            case '?':
                return EXIT_FAILURE;

            default:
                abort();
        }
    }

    char nucleotides[5] = {'A', 'C', 'G', 'T', 'N'};
    double nuc_cs[5] = {0.28, 0.49, 0.70, 0.90, 1.00};

    char qualities[64];
    double qual_cs[64];
    size_t i;
    double last_c = 0.0;
    for (i = 0; i < 64; ++i) {
        qualities[i] = '!' + i;
        last_c = qual_cs[i] = last_c + 1.0 / 64.0;
    }

    char id_chars[94];
    double id_cs[94];
    last_c = 0.0;
    for (i = 0; i < 94; ++i) {
        id_chars[i] = '!' + i;
        last_c = id_cs[i] = last_c + 1.0 / 94.0;
    }

    char* id   = malloc(id_len + 1);
    char* seq  = malloc(max_len + 1);
    char* qual = malloc(max_len + 1);
    size_t len = min_len;

    while (1) {
        if (max_len > min_len) {
            len = min_len + (size_t) (drand48() * (double) (1 + max_len - min_len));
        }

        randcat(id_chars, id_cs, 94, id, id_len);
        id[id_len] = '\0';

        randcat(nucleotides, nuc_cs, 5, seq, len);
        seq[len] = '\0';

        randcat(qualities, qual_cs, 64, qual, len);
        qual[len] = '\0';

        printf(
            "@%s\n%s\n+\n%s\n",
            id, seq, qual);
    }

    /* Yeah, right. As if we'll ever get here. */
    free(id);
    free(seq);
    free(qual);

    return EXIT_SUCCESS;
}


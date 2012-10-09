
#include "../src/parse.c"
#include "../src/common.c"

int main()
{
    seq_t* seq = seq_create();
    fastq_t* f = fastq_create(stdin);

    while (fastq_read(f, seq)) {
        fastq_print(stdout, seq);
    }

    fastq_free(f);
    seq_free(seq);

    return EXIT_SUCCESS;
}




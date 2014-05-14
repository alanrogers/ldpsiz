#include <stdio.h>
#include <getopt.h>

static int  optndx;
static int  verbose_flag = 0;

static struct option myopts[] = {
    /* {char *name, int has_arg, int *flag, int val} */
    {"verbose", no_argument, &verbose_flag, 1}, /* set verbose_flag=1 */
    {"noarg", no_argument, 0, 'n'}, /* return 'n' */
    {"file", required_argument, 0, 'f'},    /* return 'f' */
    {"maybe", optional_argument, 0, 'm'},   /* return 'm' */
    {NULL, 0, NULL, 0}          /* last element all zeroes */
};

int main(int argc, char **argv) {
    int         i;

    /*
     * According to the man page, getopt_long_only returns ':' for
     * options that require a value but lack one. In reality, the
     * behavior depends on whether there is a subsequent option on the
     * command line. If an option that requires a value is immediately
     * followed by another option (e.g. "--file --maybe"), then
     * "--maybe" is erroneously interpreted as the value of
     * "--file". On the other hand, if there is nothing following
     * "--file" on the command line, then the function prints an error
     * message and returns '?', the same value it returns for
     * unrecognized options.
     *
     * Options with optional arguments do not recognize those
     * arguments unless the option and its argument are connected by
     * '=' on the command line. In other words, this works:
     *
     *      --maybe=foo
     *
     * but this does not:
     *
     *      --maybe foo
     *
     * In the optstring argument, "vnf:m::ox", I have provided a
     * single-character alternative for each of the full-word options
     * listed in myopts. The colon after 'f' says that -f requires an
     * option. The double colon after 'm' says that -m takes an
     * optional argument.
     */
    for(;;) {
        i = getopt_long_only(argc, argv, "vnf:m::ox", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            fprintf(stderr, "error in command line arguments\n");
            break;
        case 0:
            printf("option %s\n", myopts[optndx].name);
            break;
        case 'v':
            printf("got v\n");
            break;
        case 'n':
            printf("got n\n");
            break;
        case 'f':
            printf("got f\n");
            if(optarg != NULL)
                printf("  val=%s\n", optarg);
            else
                printf("  this should not happen\n");
            break;
        case 'm':
            printf("got m\n");
            if(optarg != NULL)
                printf("  val=%s\n", optarg);
            else
                printf("  optional argument not provided for \"%s\"\n",
                       myopts[optndx].name);
            break;
        case 'x':
            printf("got x optndx=%d\n", optndx);
            printf("  optarg=%lu\n", (long unsigned) optarg);
            if(optarg != NULL)
                printf("  val=%s\n", optarg);
            break;
        case ':':
            /* According to the docs, we should get here if no
             * argument is provided to an option that requires one.
             * This doesn't seem to work.
             */
            printf("missing argument; optndx=%d\n", optndx);
            break;
        default:
            printf("?? getopt_long_only returned [%c]\n", i);
        }
    }
    if(optind < argc) {
        printf("non-flag args:");
        while(optind < argc)
            printf(" %s", argv[optind++]);
        putchar('\n');
    }
    return 0;
}

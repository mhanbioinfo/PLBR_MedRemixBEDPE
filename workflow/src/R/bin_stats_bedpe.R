# docopt
' bin_stats v1.0
Compute BAM stats by bin

Usage:
    bin_stats_bedpe.R -i INPUT -g GENOME -c CHROM -o OUT_FILE [ --bsgchr BSGCHR --winsize WINDOWSIZE ]

Options:
    -i --input INPUT        Path to cleaned BEDPE4MEDREMIX file (full path)
    -g --genome GENOME      Either the name of BSgenome (usually BSgenome.Hsapiens.UCSC.hg38
                                or BSgenome.Athaliana.TAIR.TAIR9), or the path to a folder
                                containing a custom BSgenome as a package, which will be loaded
                                using devtools::load_all()
    -c --chrom CHROM        Chromosome
    -o --outfile OUT_FILE   Path to output bin_stat file (full path)

    --bsgchr BSGCHR         If CHROM does not match the corresponding chromosome name in BSgenome,
                                provide here the actual name of the corresponding BSgenome chromosome
    --winsize WINDOWSIZE    Window size

' -> doc

## --filtered FILTERED     Path to dump reads that were filted out of the analysis
## currently not yet supported

# library
if (! interactive()) {
    suppressMessages(library(docopt))
    args <- docopt(doc, version='Bin Stats vs. 1.0')
    print(args)
} else {
    message('Running in interactive mode. Be sure to specify args manually.')
}
suppressMessages(library(arrow))
suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome))
suppressMessages(library(Rsamtools))
suppressMessages(library(tidyverse))


if (!is.null(args[['winsize']])) {
    BIN_WIDTH = args[['winsize']] %>% as.integer() 
} else {
    BIN_WIDTH = 300
}

FRAGMENT_LENGTH_LIMIT = 500

if (file.exists(paste(args[['genome']], 'DESCRIPTION', sep='/'))) {
    devtools::load_all(args[['genome']])
    bsgenome <- getBSgenome(basename(args[['genome']]))
} else {
    bsgenome <- getBSgenome(args[['genome']])
}

if (!is.null(args[['bsgchr']])) {
    stopifnot(args[['bsgchr']] %in% seqnames(bsgenome))
    bsgenome_chr = args[['bsgchr']]
} else {
    bsgenome_chr = args[['chrom']]
}

chrom_length <- length(bsgenome[[args[['chrom']]]])
bins = GRanges(
    seqnames = args[['chrom']],
    ranges = IRanges(
        start = seq(1, chrom_length, BIN_WIDTH),
        end = c(seq(BIN_WIDTH, chrom_length, BIN_WIDTH), chrom_length)
    )
)


# read in bedpe4medremix file
message(paste0("Reading in bedpe4medremix file: ", args[['input']]))

fragments_tibble_limited =
  read.table(file = args[['input']], header = F, sep = '\t') %>%
  dplyr::rename("seqnames" = "V1", "mateid" = "V2", "paired_end_reads" = "V3", "start" = "V4",
         "end" = "V5", "width" = "V6", "mean_mapq" = "V7") %>%
  dplyr::mutate(mateid = 0)


# find overlaps with bins
fragment_bin_overlaps <- findOverlaps(GRanges(fragments_tibble_limited), bins)

bin_coverage <- bind_cols(
        fragments_tibble_limited[queryHits(fragment_bin_overlaps), ],
        bins[subjectHits(fragment_bin_overlaps)] %>%
            as_tibble %>%
            rename(bin_chr = seqnames, bin_start = start, bin_end = end) %>%
            select(-width)
    ) %>%
    mutate(
        overlap_length = pmin(end, bin_end) - pmax(start, bin_start) + 1
    ) %>%
    group_by(bin_chr, bin_start, bin_end) %>%
    summarise(
        n_fragments = n(),
        coverage_bp = sum(overlap_length),
        mean_fragment_length = mean(width),
        mean_fragment_mapq = mean(mean_mapq)
    ) %>%
    ungroup() %>%
    right_join(
        bins %>%
            as_tibble %>%
            select(seqnames, start, end),
        by = c('bin_chr' = 'seqnames', 'bin_start' = 'start', 'bin_end' = 'end')
    ) %>%
    replace_na(list(
        n_fragments = 0,
        coverage_bp = 0
    )) %>%
    arrange(bin_chr, bin_start) %>%
    mutate(mean_coverage = coverage_bp / BIN_WIDTH) %>%
    mutate(
        seq = getSeq(
            bsgenome,
            names = rep(bsgenome_chr, n()),
            start = bin_start,
            end = bin_end
        ) %>% as.character
    ) %>%
    mutate(
        known_bps = str_count(seq, '[TCGA]'),
        gc_content = str_count(seq, '[GC]') / ifelse(known_bps == 0.0, 1.0, known_bps),
        cpg_count = str_count(seq, 'CG'),
    ) %>%
    select(-known_bps, -seq)
# bin_coverage

bin_coverage %>%
    write_tsv(
        args[['outfile']],
        append=TRUE,
        col_names = TRUE
    )

if (!is.null(args[['filtered']])) {
    filtered_reads %>% write_tsv(args[['filtered']])
}



# EOF

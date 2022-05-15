# docopt
' bin_stats_bam v1.0
Compute BAM stats by bin

Usage:
    bin_stats_bam.R -i INPUT -g GENOME -c CHROM -o OUTPUT [ --bsgchr BSGCHR --filtered FILTERED --winsize WINDOWSIZE ]

Options:
    -i --input INPUT        Path to input BAM or BEDPE file
    -g --genome GENOME      Either the name of BSgenome (usually BSgenome.Hsapiens.UCSC.hg38
                                or BSgenome.Athaliana.TAIR.TAIR9), or the path to a folder
                                containing a custom BSgenome as a package, which will be loaded
                                using devtools::load_all()
    -c --chrom CHROM        Chromosome
    -o --output OUTPUT      Output path

    --bsgchr BSGCHR         If CHROM does not match the corresponding chromosome name in BSgenome,
                                provide here the actual name of the corresponding BSgenome chromosome
    --filtered FILTERED     Path to dump reads that were filted out of the analysis
    --winsize WINDOWSIZE    Window size

' -> doc


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


## process BAM
message(paste0("Processing: ", args[['input']], "."))
bam_file <- BamFile(args[['input']], asMates=TRUE)
bam_reads <- scanBam(
  bam_file,
  param = ScanBamParam(
    which = GRanges(sprintf('%s:1-%s', args[['chrom']], chrom_length)),
    what = scanBamWhat()
  )
)

message("Getting bam_reads_tibble.")
bam_reads_tibble <- tibble(
  seqnames = args[['chrom']],
  start = bam_reads[[1]]$pos,
  width = bam_reads[[1]]$qwidth,
  mate_status = bam_reads[[1]]$mate_status,
  mateid = bam_reads[[1]]$groupid,
  strand = bam_reads[[1]]$strand,
  mapq = bam_reads[[1]]$mapq,
  flag = bam_reads[[1]]$flag
) %>%
  mutate(end = start + width) %>%
  dplyr::filter(!(str_detect(flag, "2[0-9][0-9][0-9]") & mate_status == "mated"))
## { deal with extra reads where flag says SUPPLEMENTARY but still has mate_status mated }


## deal with these freq3 alignments
## leave freq4 etc as ambiguous and removed later as usual

bam_reads_tibble.mateid.table =
  bam_reads_tibble$mateid %>% table() %>% as.data.frame()

bam_reads_tibble.mateid.table.freq3 =
  bam_reads_tibble.mateid.table %>%
  dplyr::filter(Freq == 3) %>%
  dplyr::rename("mateid" = ".") %>%
  dplyr::mutate(mateid = as.character(mateid))

bam_reads_tibble.no_mateid_freq3 =
  bam_reads_tibble %>%
  dplyr::filter(!mateid %in% bam_reads_tibble.mateid.table.freq3$mateid)

## distinct just for 'seqnames, start, mateid, strand', cuz 'end' sometimes differ...
## want to collapse these freq3 reads into freq2 so matches .bedpe4medremix
mateid_freq3_fixed =
  bam_reads_tibble %>%
  dplyr::filter(mateid %in% bam_reads_tibble.mateid.table.freq3$mateid) %>%
  dplyr::arrange(mateid) %>%
  dplyr::filter(!str_detect(flag, "2[0-9][0-9][0-9]")) %>%
  dplyr::distinct(seqnames, start, mateid, strand, .keep_all = TRUE) %>%
  dplyr::mutate(mate_status = "mated")

bam_reads_tibble.mateid_freq3_fixed =
  dplyr::bind_rows(bam_reads_tibble.no_mateid_freq3,
                   mateid_freq3_fixed) %>%
  dplyr::select(-flag)


# filtering
message("Getting fragments_tibble.")
fragments_tibble <-
  bam_reads_tibble.mateid_freq3_fixed %>%
  filter(mapq >= 20) %>%                      ## ANY of the mate in pair if mapq < 20, exclude
  group_by(seqnames, mateid) %>%
  filter(n() == 2, '-' %in% strand, '+' %in% strand) %>%    ## eliminate any reads with >2 mateids
  filter(mate_status == 'mated') %>%
  summarise(
    paired_end_reads = sprintf(
      '+:%s-%s,-:%s-%s',
      start[strand=='+'],
      end[strand=='+'],
      start[strand=='-'],
      end[strand=='-']
    ),
    start = start[strand == '+'],
    end = end[strand == '-'],
    width = end-start,
    mean_mapq = mean(mapq)
  ) %>% ungroup

fragments_tibble_limited <- fragments_tibble %>% filter(width < FRAGMENT_LENGTH_LIMIT, width > 0)

message(paste0("Finished processing ", args[['chrom']], " fragment tibble."))

filtered_reads <- bam_reads_tibble %>% filter(! mateid %in% fragments_tibble_limited$mateid)

message(
    sprintf('Filtered out %s fragments with insufficient quality, length <= 0, or length > %s',
    nrow(fragments_tibble) - nrow(fragments_tibble_limited),
    FRAGMENT_LENGTH_LIMIT
    )
)


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

bin_coverage %>%
    write_tsv(
        args[['output']],
        append=TRUE,
        col_names = TRUE
    )

if (!is.null(args[['filtered']])) {
    filtered_reads %>% write_tsv(args[['filtered']])
}



# EOF

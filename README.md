# filter-clipped #

[Rust learning project]

Remove alignments with high number of clipped base.

Sometimes aligners have very loose scoring methods and write alignments with high abundant of soft/hard-clipped base into alignment BAM files. This program is for filtering these reads out by gating the number of clipped bases in relative to the read sequence length

## Installation


```
$ git clone https://github.com/wckdouglas/filter-clipped
$ cd filter-clipped
$ cargo install --path .  # if compilation error, try CC=/usr/bin/gcc cargo install --path .
$ filter-clipped --help
filter-clipped 0.1.0
Remove alignments with high number of clipped base. Sometimes aligner has very loose scoring methods
and write alignments with high abundant of soft/hard-clipped base into alignment BAM files. This
program is for filtering these reads out by gating the number of clipped bases in relative to the
read sequence length

USAGE:
    filter-clipped [OPTIONS] --in-bam <IN_BAM>

OPTIONS:
    -b, --both-end <BOTH_END>        maximum fraction of total bases on the sequence being clipped
                                     [default: 0.1]
    -h, --help                       Print help information
    -i, --in-bam <IN_BAM>            input bam file path  ("-" for stdin)
        --inverse                    keeping the failed ones (high-clipped-fraction alignments)
    -l, --left-side <LEFT_SIDE>      maximum fraction of bases on the sequence being clipped from
                                     the left side (5' end) [default: 0.1]
    -o, --out-bam <OUT_BAM>          output bam file path ("-" for stdout) [default: -]
    -r, --right-side <RIGHT_SIDE>    maximum fraction of bases on the sequence being clipped from
                                     the right side (3' end) [default: 0.1]
    -V, --version                    Print version information
```


## Test 
```
cargo test
```

## Docker 
```
docker pull ghcr.io/wckdouglas/filter-clipped:main
docker run --env RUST_LOG=info  ghcr.io/wckdouglas/filter-clipped:main --in-bam test/data/test.sam | samtools view
```


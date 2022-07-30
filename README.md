# filter-clipped #

Remove alignments with high number of clipped base.

Sometimes aligners have very loose scoring methods and write alignments with high abundant of soft/hard-clipped base into alignment BAM files. This program is for filtering these reads out by gating the number of clipped bases in relative to the read sequence length

## Installation


```
git clone https://github.com/wckdouglas/filter-clipped
cd filter-clipped
cargo install # if compilation error, try CC=/usr/bin/gcc cargo install
filter-clipped --help
```


## Test 
```
cargo test
```
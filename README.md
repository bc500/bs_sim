# bs_sim
## Bisulfite treatment simulator for DNA sequences

Bisulfite treatment is a protocol used in sequencing to detect the methylation state of cytosines in a DNA sequence. It replaces unmethylated cytosines by thymines and, after sequencing and alignment, this reveals the methylation state of a genomic position.

The bisulfite treatment simulator takes paired-end sequencing reads and, using a methylation probability profile, simulates the sample preparation using a treatment of bisulfite. This generates a modified DNA sequence with converted cytosines. The conversion efficiency rates can be adjusted using the parameters λ (unmethylated cytosine conversion rate) and τ (methylated cytosine conversion rate).

## Installation
* Download the code 
```
git clone https://github.com/bc500/bs_sim.git
```
* Compile
```
cd bs_sim
./configure && make
```
This will generate the binary ```src/bs_sim```.

* Install

The file can be either copied into your user bin directory
```
cp src/bs_sim ~/bin/
```
or installed system-wide
```
sudo make install
```


## Usage
```
Usage:
  bs_sim --help
  bs_sim --version
  bs_sim -i INPUT.BAM -o OUTPUT.BAM [--lambda=<rate>] [--tau=<rate>] 
          [--methylation-profile=<PROFILE.TBX>]
          [--cpg-islands=<CPGIslands.TBX> --methylation-cpg=<rate>]
          [--seed=<integer>]

Options:
  -h --help           Show this screen.
  -v --version        Show version.
  -i --input=FILE     Input SAM/BAM file.
  -o --output=FILE    Output SAM/BAM file.
  -l --lambda=<rate>  Probability of conversion from cytosine
                      to thymine [default: 0.99].
  -t --tau=<rate>     Probability of conversion from methylated cytosine
                      to thymine [default: 0.05].
  -p --methylation-profile=<PROFILE.TBX>
                      Tabix-format file with methylation probability for 
                      each genomic position.
  -c --cpg-islands=<CPGIslands.TBX>
                      Tabix-format file with CpG Island regions.
  -m --methylation-cpg=<rate>
                      Probability of cytosine being methylated on a CpG Island [default: 0.90].
  -s --seed=<integer> Seed for random number generation,
                      to make experiments repeatable.
```
## Examples
```

./bs_sim -i sample.bam -o sample_bs.bam -p methylation.tbx -c CpGs.tbx

./bs_sim -i sample.bam -o sample_bs.bam -p methylation.tbx -c CpGs.tbx -l 0.98 -t 0.03 -m 0.8

```

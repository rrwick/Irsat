# Irsat: Iterative Read Subset Assembly Tool

Irsat is a Python script that conducts a focused, iterative *de novo* assembly.

It works by first mapping reads to a target sequence, then conducting *de novo* assembly on the reads that mapped as well as their read pairs. The process then repeats, now with the assembled contigs added to the mapping reference, for as many iterations as necessary.

Each iteration generates the following:
* Read subset files (FASTQ)
* Assembled contigs (FASTA)
* Assembly graph

Ideally, the iteration will grow with each iteration. The resulting assembly graph can then be analysed using Bandage.

## Requirements

* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Samtools](https://samtools.github.io/)
* [Bedtools](http://bedtools.readthedocs.io/)
* An assembler, e.g. [SPAdes](http://cab.spbu.ru/software/spades/)

## Installation

No compilation or installation is required - just download/clone and run Irsat.py.

## Contributing

New contributors are welcome!  If you're interested or have ideas, please contact me (Ryan) at rrwick@gmail.com.

## History

Version 0.1.0 – initial release on GitHub

Version 0.1.1 – minor updates to handle current versions of Samtools and SPAdes.

## License

GNU General Public License, version 3

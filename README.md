# Introduction

CountMap is a C/C++ implementation of a short read mapper. For the alignment phase it uses the KSW2 algorithm implemented by [Heng Li][hl]. The implementation was made as a part of my Final BSc thesis.

# Dependencies

1. gcc 4.8+
2. cmake 3.2+

# Installation

In order to build CountMap from source run the following commands:

```bash
git clone --recursive --single-branch -b spavlic https://github.com/lbcb-edu/BSc-thesis-18-19.git countmap
cd countmap
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

After running the commands, the executable file will be created in `build/bin`.

# Usage

CountMap outputs mappings in the [SAM format][sam].

For mapping single-end reads:
```bash 
countmap reference.fa query.fq > single-end_mappings.sam
```

For mapping paired-end reads:
```bash
countmap -p reference.fa query1.fq query2.fq > paired-end_mappings.sam
```

For mapping paired-end reads as single-end:
```bash
countmap reference.fa query1.fq query2.fq > paired-end_single_mappings.sam
```

For advanced options check the `help` page with:
```bash
countmap --help
```

# Input data test examples

Escherichia coli str. K-12 substr. MG1655
  - [reference genome][ref]
  - [Illumina paired-end reads][reads]

# Final BSc thesis (computer science - 2018/2019)

Final BSc thesis is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the sixth semester of the undergraduate study. The main focus is to apply knowledge and skills obtained from Software Design Project course to recreate or improve existing methods which are widely used in bioinformatics. Under the supervision of prof. Mile Šikić, students will implement one such algorithm, thoroughly test it on simulated and real data, and formally encapsulate the whole process by writing and defending a thesis. Each student will have access to a private branch of this repository, on which this README will be updated with the specific task.

## Disclaimer

Laboratory for Bioinformatics and Computational Biology cannot be held responsible for any copyright infringement caused by actions of students contributing to any of its repositories. Any case of copyright infringement will be promptly removed from the affected repositories and reported to appropriate faculty organs.

[hl]: https://github.com/lh3/ksw2
[sam]: https://samtools.github.io/hts-specs/SAMv1.pdf
[ref]: https://www.ncbi.nlm.nih.gov/genome/167
[reads]: http://www.ebi.ac.uk/ena/data/view/ERA000206&display=html
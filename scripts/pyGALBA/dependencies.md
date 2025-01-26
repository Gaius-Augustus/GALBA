## 1. Python 3 + Required Python Packages

- **Python 3** (the script is `#!/usr/bin/env python3`)
- **BioPython** (used for FASTA/FASTQ parsing: `from Bio import SeqIO`)
- **pandas** (used for TSV reading and DataFrame merging: `import pandas as pd`)
- **PyYAML** (used for reading config: `import yaml`)
- **Pygustus** (used for AUGUSTUS gene prediction: `import pygustus`)
- Built‐in libraries:
  - `argparse`, `subprocess`, `sys`, `os`, `re`, `math`, `shutil`, etc.
  - All are part of the Python standard library (no extra installation needed besides Python itself).

---

## 2. External Command‐Line Tools

The script calls the following binaries via `subprocess.run(...)`. You must install each one and ensure the script can find them (either via PATH or by specifying a directory).

1. **HISAT2**  
   - Commands used:  
     - `hisat2-build` (for indexing)  
     - `hisat2` (for short-read mapping)

2. **Minimap2**  
   - Command used:
     - `minimap2` (for mapping long reads, e.g., Iso-Seq)

3. **samtools**  
   - Commands used:
     - `samtools sort`
     - `samtools merge`

4. **StringTie**  
   - Command used:
     - `stringtie` (for assembling RNA-seq or Iso-Seq reads into transcripts)

5. **TransDecoder**  
   - Commands used:
     - `TransDecoder.LongOrfs`
     - `TransDecoder.Predict`
     - Utility scripts:  
       - `gtf_genome_to_cdna_fasta.pl`
       - `gtf_to_alignment_gff3.pl`
       - `cdna_alignment_orf_to_genome_orf.pl`

6. **DIAMOND**  
   - Commands used:
     - `diamond makedb`
     - `diamond blastp`

7. **bedtools**  
   - Command used:
     - `bedtools coverage` (for checking coverage overlaps)

8. **miniprot**  
   - Command used:
     - `miniprot` (for protein-to-genome alignments, in addition to or instead of `minimap2`)

9. **miniprot_boundary_scorer**  
   - Command used:
     - `miniprot_boundary_scorer` (for scoring protein alignment boundaries with a chosen matrix)

10. **TSEBRA**
      - Command used:
         - `tsebra` (for merging gene sets)

11. **AUGTUSTUS**
        - Command used:
             - `augustus` (for gene prediction)
             - `etraining` (for training AUGUSTUS)
                - `optimize_augustus.pl` (for optimizing AUGUSTUS)
                - `new_species.pl` (for creating a new AUGUSTUS species)

---
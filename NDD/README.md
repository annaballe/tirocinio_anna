Toolkit for benchmarking, merging, and annotating Structural Variants

## 💻 Installation
Truvari uses Python 3.6+ and can be installed with pip:
```
  python3 -m pip install Truvari 
```
For details and more installation options, see [Installation](https://github.com/acenglish/truvari/wiki/Installation) on the wiki.

## ⏩ Quick Start

Each sub-command contains help documentation. Start with `truvari -h` to see available commands.

The current most common Truvari use case is for structural variation benchmarking:
```
  truvari bench -b base.vcf.gz -c comp.vcf.gz -f reference.fa -o output_dir/
```

Find more matches by harmonizing phased varians using refine:
```
   truvari refine output_dir/
```

Use Truvari's comparison engine to consolidate redundant variants in a merged multi-sample VCF:
```
    bcftools merge -m none sampleA.vcf.gz sampleB.vcf.gz | bgzip > merge.vcf.gz
    tabix merge.vcf.gz
    truvari collapse -i merge.vcf.gz -o truvari_merge.vcf
```

# ig_cluster_splitter
ig_cluster_splitter is a tool created as one more module of igquast tool from IgRepertoireConstructor 
(https://github.com/yana-safonova/ig_repertoire_constructor). 

ig_cluster_splitter is used to divide clusters of read's collected from
.rca files and create's new .rca file of modifyed cluster's.

## Code Examples

```
python ig_cluster_splitter.py -s input_reads.fa -r igrec_output.rcm -o dir/ -f filename
```

## Motivation

Original IgReC tool suffer from overcorrection problem - mutations in antibodies often are threated as amplification errors 
and IgReC tool trying to fix them. 

In this project was made an attempt to construct ML model to reduce this effect. 

Slides about the topic and some information about results of the project are [here](https://docs.google.com/presentation/d/17QOT-wQAiNQqK-YnIj34VaYi4yobI6B1sTKvbnzlsHY/edit#slide=id.g238dd5aaef_0_7) or in [this](https://drive.google.com/open?id=0B8-bo9EAqeQgaW8wZFlCeHE3Smc) small report.

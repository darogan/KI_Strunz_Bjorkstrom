# Continuous human uterine NK cell differentiation in response to endometrial regeneration and pregnancy

**Benedikt Strunz<sup>1,*</sup>, Jonna Bister<sup>1</sup>, Russell S. Hamilton<sup>2</sup>, Hanna Jönsson<sup>1</sup>, Ylva Crona-Guterstam<sup>1,3</sup>, Egle Kvedaraite<sup>1,3</sup>, Iva Filipovic<sup>1</sup>, Natalie Sleiers<sup>1</sup>, Bogdan Dumitrescu<sup>4</sup>, Danielle Friberg<sup>5</sup>, Mats Brännström<sup>6</sup>, Antonio Lentini<sup>7</sup> , Björn Reinius<sup>7</sup>, Martin Cornillet<sup>1</sup>, Tim Willinger<sup>1</sup>, Sebastian Gidlöf<sup>3</sup>, Martin A. Ivarsson<sup>1</sup> & Niklas K. Björkström<sup>1,*</sup>**

<sup>1</sup> Center for Infectious Medicine, Department of Medicine Huddinge, Karolinska Institutet, Karolinska University Hospital, Stockholm, Sweden <br>
<sup>2</sup> Centre for Trophoblast Research, University of Cambridge, Cambridge, UK <br>
<sup>3</sup> Department of Women’s and Children’s Health, Karolinska Institutet, Stockholm, Sweden <br>
<sup>4</sup> Department of Obstetrics and Gynecology, Mälarsjukhuset, Eskilstuna, Sweden <br>
<sup>5</sup> Department of Surgical Sciences, Uppsala University, Uppsala, Sweden <br>
<sup>6</sup> Department of Obstetrics and Gynecology, University of Gothenburg, Gothenburg, Sweden <br>
<sup>7</sup> Department of Medical Biochemistry and Biophysics, Karolinska Institutet, Stockholm, Sweden <br>
<sup>*</sup> Corresponding authors: benedikt.strunz@ki.se, niklas.bjorkstrom@ki.se


### Publication ###

Strunz, B., Bister, J., Hamilton, R.S., Jönsson, H., Crona-Guterstam, Y., Kvedaraite, E.,  Filipovic, I., Sleiers, N., Dumitrescu, B., Friberg, D., Brännström, M., Lentini, A., Reinius, B., Cornillet, M., Willinger, T., Gidlöf, S., Ivarsson, M.A., & Björkström, N.K. [[<s>JOURNAL</s>]](https://) [[<s>DOI</s>]](https://doi.org/) [[<s>bioRxiv</s>]](https://doi.org/10.1101/)

### Abstract ###

On publication

### Bulk RNA-Seq ###

#### Sample Table ####

| Sample Name	| Group      |
| ----------- | ---------- |
| 1_HU12_P11  | CD39+KIR+  |
| 2_HU12_P9	  | CD39-KIR-  |
| 3_HU12_P10  | CD39-KIR+  |
| 5_HU13_P11  | CD39+KIR+  |
| 6_HU13_P9	  | CD39-KIR-  |
| 7_HU13_P10  | CD39-KIR+  |
| 9_Hu05_P11	| CD39+KIR+  |
| 10_Hu05_P9	| CD39-KIR-  |
| 11_Hu05_P10	| CD39-KIR+  |
| 14_Hu7_P12	| CD39+KIR+  |
| 15_Hu7_P13	| CD39-KIR+  |
| 16_Hu7_P9	  | CD39-KIR-  |


#### Data Processing ####

Raw sequencing files are run through quality control using FastQC (v0.11.5) and fastq_screen (v0.9.3). Low quality and adapter sequencing are trimmed with Trim Galore! (v0.6.4). Trimmed reads are aligned to the reference genome (GRCh38, ensEMBL) using STAR (v020201). Alignments are assessed using qualimap (v2.2) and featureCounts (v 1.5.0-p2). Gene quantification is performed with featureCounts (v 1.5.0-p2). Differential gene expression is performed with DESeq2 package (v1.22.2, R v3.5.3), including principle component analysis (PCA) to assess sample clustering, and multiple testing correction to produce false discovery rates. Finally all metrics from the RNA-Seq pipelines are summarised and reports produced using MultiQC (0.9.dev0).

> Additional Processing Steps

CV2

> Resources Used

Resource       | URL
-------------- | --------------
GRCh38         | [Link](https://www.ensembl.org/Homo_sapiens/Info/Index)
FastQC         | [Link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
Trim_galore    | [Link](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
STAR           | [DOI](https://doi.org/10.1093/bioinformatics/bts635)
HTSeq-counts   | [DOI](http://dx.doi.org/10.1093/bioinformatics/btu638)
Feature_counts | [DOI](http://dx.doi.org/10.1093/bioinformatics/btt656)
Qualimap       | [DOI](https://doi.org/10.1093/bioinformatics/bts503)
RSeQC          | [DOI](http://doi.org/10.1093/bioinformatics/bts356)
ClusterFlow    | [DOI](http://dx.doi.org/10.12688/f1000research.10335.2)
MultiQC        | [DOI](http://dx.doi.org/10.1093/bioinformatics/btw354)

##### Bulk RNA-Seq Pipeline #####

Pipeline run using [ClusterFlow](http://clusterflow.io)

    #fastqc
    #fastq_screen
    #trim_galore
          #fastqc
          #star
            #qualimap_rnaseq
            #rseqc_infer_experiment
            #featureCounts
            #htseq_counts

FeatureCount gene count files are available in the [FeatureCount](FeatureCount/) directory

#### Script to reproduce paper figures (Bulk RNA-Seq only) ####

Figure        | Description
------------- | --------------
Figure 2G     | PCA


### Session Information ###
Details for the R version and packages used to create all figure



### Links ###

Description   | URL
------------- | ----------
Publication   | [[<s>JOURNAL</s>]](https://) [[<s>DOI</s>]](https://doi.org/) [[<s>bioRxiv</s>]](https://doi.org/10.1101/)
Raw Data      | ArrayExpress EMBL-EBI [<s>E-MTAB-????</s>](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-????)
Björkström Group | [Björkström group website](https://ki.se/en/medh/niklas-bjorkstrom-group)
CTR Bioinformatics | [CTR-BFX](https://www.trophoblast.cam.ac.uk/Resources/BioInformatics)

### Contact ###

Contact Russell S. Hamilton (rsh46 -at- cam.ac.uk) for bioinformatics related queries

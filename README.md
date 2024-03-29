# Reciprocal Best Structure Hits
## Background
In this work, we are using AlphaFold structure models [[1]](#1) to find the closest homologues proteins between <i>Homo sapiens</i> and <i>D. melanogaster</i>, <i>C. elegans</i>, <i>S. cerevisiae</i> and <i>S. pombe</i> as well as between <i>S. cerevisiae</i> and <i>S. pombe</i>. We are using the strcuture aligner Foldseek [[2]](#2) to run all against all and search for the best scoring hit in both directions to detect the Reciprocal Best Structure Hits (RBSH). We compare the results to protein pairs detected by their sequence similarity as Reciprocal Best Hits (RBH) and verify the results using the PANTHER family classification files [[3]](#3). Our work is described in [[4]](#4)</br>

## Dataset
The detected RBSH and RBH can be found under [Results_file](https://github.com/VivianMonzon/Reciprocal_Best_Structure_Hits/tree/main/Results_files). 

## References
<a id="1">[1]</a>
Jumper J. et al. Highly accurate protein structure prediction with AlphaFold.
Nature, doi:[10.1038/s41586-021-03819-2](https://www.nature.com/articles/s41586-021-03819-2) (2021)

<a id="2">[2]</a>
Van Kempen M. et al. Foldseek: fast and accurate protein structure search.
bioRxiv, doi:[10.1101/2022.02.07.479398](https://www.biorxiv.org/content/10.1101/2022.02.07.479398v4) (2022)

<a id="3">[3]</a>
Mi H. et al. PANTHER version 16: a revised family classification, tree-based classification tool, enhancer regions and extensive API.
Nucleic Acids Res., doi:[10.1093/nar/gkaa1106](https://academic.oup.com/nar/article/49/D1/D394/6027812?login=true) (2021)

<a id="4">[4]</a>
Monzon V., Paysan-Lafosse T., Wood V., Bateman A. Reciprocal best structure hits: using AlphaFold models to discover distant homologues.
Bioinformatics Advances, doi:[10.1093/bioadv/vbac072](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac072/6749558) (2022)
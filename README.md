# Digenic

Instructions on how to use Rscripts are available on the Rscripts themselves. 
Scripts were developped to be used under double loop parallelisation using internal R packages. 
However it is recommended to follow the advice of manual parallelisation as stated in the scripts. 
Detail information is expected to be available in the near future as well as updated scripts.

# General information

Scripts were developed to perform digenic analyses. The script example is shown for gene x gene interaction analyses. However gene x variant is as straightforward as gene x gene interaction tests using these scripts. 
For gene x gene interaction tests, the general idea is to aggregate variant information at the gene level and assess unexpected high correlation between carriers of variants of two different genes. The test in conducted under a case-only framework which provides the advantage of rely exclusively on absence of population stratification and LD to assess true interactions. 

DiversiTools
============

Test set1: avian influenza H7N1
-------------------------------

These samples come from the NGS study of avian H7N1 outbreak in 1999
<ul>
<li>S20 LPAI (4829-99)</li>
<li>S18 LPAI (4295-99)</li>
<li>S3 HPAI (4827-99)</li>
</ul>

**Reference:**

Monne I, Fusaro A, Nelson MI, Bonfanti L, Mulatti P, Hughes J, Murcia PR, Schivo A, Valastro V, Moreno A, Holmes EC, Cattoli G.
Emergence of a highly pathogenic avian influenza virus from a low-pathogenic progenitor. J Virol. 2014 Apr;88(8):4375-88. 
doi: 10.1128/JVI.03181-13
PMID: 24501401 

**Example:**

    
    perl diversiutils.pl -bam test/S3_refHPAI_cons_stampy.bam -orf test/Coding.regions.HPAI.txt -ref test/refHPAI_cons.fa
    

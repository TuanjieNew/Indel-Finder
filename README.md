# Indel-Finder  
A program used to find indels from next-genearation sequencing data on target gene edited by CRISPR/Cas9 and then calculate the on-target efficiency. The indel information will be presented as figures plotted by matplotlib, a module of python. You can also use other software, like OriginLab, to plot indel data by yourself. 

## Environment configuration:  
<code>Indel-Finder</code> uses [Nucleotide-Nucleotide BLAST 2.4.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) to find reads position and do local alignment and then processes the blast result to get indel and calculate the on-target efficiency.  


## Dependency:  
<code>Indel-Finder</code> use numpy module and matplotlib module.  
User can install them by using pip:

  > pip install numpy  
  > pip install matplotlib  

## Get Indel-Finder:  
> git clone https://github.com/TuanjieNew/Indel-Finder.git  

## Usage: 
>python indel_finder.py -r ./sample/crr5ref.fa -t ./sample/target1.fa -1 ./sample/1_S1_L001_R1_001.fastq -2 ./sample/1_S1_L001_R2_001.fastq -o ./    
 
## Full options:  
### *Common options*  
<code> -h &emsp;&emsp;show this help message and exit.</code>  
### *File options*  
<code> -r &emsp;&emsp;file name of reference file, required. The file is fasta format.</code>    
<code> -t &emsp;&emsp;file name of target sequence, required. The file is fasta format.</code>   
<code> -1 &emsp;&emsp;file name of read1, required. The file is fastq or fastq.gz format.</code>    
<code> -2 &emsp;&emsp;file name of read2, required. The file is fastq or fastq.gz format.</code>    
<code> -o &emsp;&emsp;outpout directories. Default is ./ . The blast results are in cas_temp, indel data in ./*_result, figures in ./*_figures.</code>  

## Mosaicism Index  

Mosacism index is from [Gini–Simpson index](https://en.wikipedia.org/wiki/Diversity_index) used to evaluate the degree of chimera.  
A mosacism index close to 1 means a high degree of chimera.  
It can be  expressed as a transformation of true diversity of order 2:  
![Gini-Simpson index](https://wikimedia.org/api/rest_v1/media/math/render/svg/cfe79cc21d7d7f882b22f2ef6660ba2461640246)


## License  
<code>Indel-Finder</code> is [MIT-licensed](https://github.com/TuanjieNew/Indel-Finder/blob/master/LICENSE).  

# Indel-Finder  
A program used to find indels from next-genearation sequencing data on target gene edited by CRISPR/Cas9 and then calculate the on-target efficiency.  

## Environment configuration:  
<code>Indel-Finder</code> uses [Nucleotide-Nucleotide BLAST 2.4.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) to find reads position and do local alignment and then processes the blast result to get indel and calculate the on-target efficiency.  


## Depedency:  
<code>Indel-Finder</code> use numpy module and matplotlib module.  
User can install them by using pip:

  > pip install numpy  
  > pip install matplotlib  

## Get Indel-Finder:  
> git clone https://github.com/TuanjieNew/Indel-Finder.git  

## Usage: 
>python bindel.py -r ./sample/crr5ref.fa -t ./sample/target1.fa -1 ./sample/1_S1_L001_R1_001.fastq -2 ./sample/1_S1_L001_R2_001.fastq -o ./    
 
## Full options:  
### *Common options*  
<code> -h &emsp;&emsp;show this helop message and exit.</code>  
### *File options*  
<code> -r &emsp;&emsp;file name of reference file, required. The file is fasta format.</code>    
<code> -t &emsp;&emsp;file name of target sequence, required. The file is fasta format.</code>   
<code> -1 &emsp;&emsp;file name of read1, required. The file is fastq or fastq.gz format.</code>    
<code> -2 &emsp;&emsp;file name of read2, required. The file is fastq or fastq.gz format.</code>    
<code> -o &emsp;&emsp;outpout directories. Default is ./</code>  


## License  
<code>Indel-Finder</code> is [MIT-licensed](https://github.com/TuanjieNew/Indel-Finder/blob/master/LICENSE).  

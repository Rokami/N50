# N50
 usage:perl AsmN50.pl <scaffold.fasta> [repeat.gff3]    
>if only supply scaffold.fasta,exec will calculate ATGCN percentage and N50;    
>if reapeat.gff3 supplied,exec will calculate percentage of Repeat percentage and output    
>overlap of each class of repeat to log.txt.   
>and transfer Repeat sequence to lowchar,output to <scaffold.fasta>.transfer.fasta;   
            
`scaffold.fasta`:Assembled fasta file  
`repeat.gff3`:generate by RepeatMasker

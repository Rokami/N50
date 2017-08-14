#!/usr/bin/perl
use Bio::SeqIO;
if(!$ARGV[0])
{
    print <<END;
    usage:perl AsmN50.pl <scaffold.fasta> [repeat.gff3]
          if only supply scaffold.fasta,exec will calcuate ATGCN percentage and N50;
          if reapeat.gff3 supplied,exec will calcuate percentage of Repeat percentage and output 
          overlap of each class of repeat to log.txt.
          and transfer Repeat sequence to lowchar,output to <scaffold.fasta>.transfer.fa;
END
exit(256);
}
$in=Bio::SeqIO->new(-file=>$ARGV[0],-format=>'fasta');
while($seq=$in->next_seq)
{
            $fasta{$seq->display_name}=$seq->seq;
                (@base)=split //,$seq->seq;
                    push @n,$seq->length;
                    foreach $key(@base)
                    {
                        $hash{$key}++;
                    }

}
grep{$sum+=$_}values %hash;
print "A:$hash{A}\t",sprintf("%.2f%%",100*$hash{A}/$sum),"\n";
print "G:$hash{G}\t",sprintf("%.2f%%",100*$hash{G}/$sum),"\n";
print "C:$hash{C}\t",sprintf("%.2f%%",100*$hash{C}/$sum),"\n";
print "T:$hash{T}\t",sprintf("%.2f%%",100*$hash{T}/$sum),"\n";
print "N:$hash{N}\t",sprintf("%.2f%%",100*$hash{N}/$sum),"\n";
undef %hash;
@n=sort{$b<=>$a}@n;
$flag=-1;
$count=0;
foreach(@n)
{
    $count+=$_;
    if($count>$sum/2)
    {
        print "N50:$n[$flag]\n";
        last;
    }
    $flag++;
}
$count=0;
exit(0) if(!$ARGV[1]);
open FN,$ARGV[1] or die $!;
while(<FN>)
{
    @line=split /\t/;
    $t=substr($fasta{$line[0]},($line[3]-1),($line[4]-$line[3]+1));
    substr($fasta{$line[0]},$line[3]-1,$line[4]-$line[3]+1,lc $t);
    $count+=($line[4]-$line[3]+1);
    $temp=$line[0] if(!$temp);
    @info=split /;/,$line[-1];
    $info[2]=~/Class=(.*)/;
    if($line[0] ne $temp)
    {
        @log=sort{${split ",",$a}[1] cmp ${split ",",$b}[1]}@log;
        foreach(@log)
        {
            @line_1=split ",";
            $type=$line_1[1]if(!$type);
            if($type ne $line_1[1])
            {
                @pos=sort{$a<=>$b}@pos;
                @p=@pos;
                while(@pos)
                {
                    $former=shift @pos;
                    @head=split ",",$former;
                    foreach $back(@pos)
                    {
                       @back=split ",",$back;
                       if($head[1]>=$back[0])
                       {
                           if($head[1]>=$back[1])
                           {
                               $hash{$back[0].",".$back[1]}=0;
                           }
                           else
                           {
                               $hash{$back[0].",".$head[1]}=0;
                           }
                       }
                    }
                }
                foreach$key(keys %hash)
                {
                    ($start,$end)=split ",",$key;
                    $length=$end-$start+1;
                    foreach $t(@p)
                    {
                        @pair=split ",",$t;
                        if($start>=$pair[0] && $end<=$pair[1])
                        {
                            $hash{$key}++;
                        }
                    }
                    push @out,"$temp\t$type\t$key\t$length\t$hash{$key}\n";
                }
                undef %hash;
                undef @pos;
            }
            push @pos,"$line_1[2],$line_1[3]";
            $type=$line_1[1];
        }
        undef @log;
    }
    push @log,"$line[0],$1,$line[3],$line[4]";
    $temp=$line[0];
}
grep{$m+=((split /\t/)[-1]-1)*(split /\t/)[-2]}@out;
print "Total TE length:",$count-$m,"\n";
print "TE Percentage of Scaffold:",sprintf("%.2f%%",100*($count-$m)/$sum),"\n";
@out=sort{(split /\t/,$b)[-1] <=> (split /\t/,$a)[-1]}@out;
open OUT,">log.txt"or die $!;
$i=(split /\t/,$out[0])[-1];
while($i>=2)
{
print OUT sort{ (split /\t/,$b)[3]<=>(split /\t/,$a)[3]}grep{(split /\t/,$_)[-1]==$i}@out;
$i--;
}
close IN;
open IN,">$ARGV[0].transfer.fasta"or die $!;
foreach $key (sort{$a cmp $b}keys %fasta)
{
    print  IN ">$key\n";
    print  IN "$fasta{$key}\n";
}
close IN;
close FN;
close OUT;

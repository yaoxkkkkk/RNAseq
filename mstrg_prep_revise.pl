#!/bin/env perl
#Usage: mstrg_prep.pl merged.gtf > merged_prep.gtf
use strict;
my %g; # gene_id => \%ref_gene_ids (or gene_names)
my @prep; # array of [line, original_id]
while (<>) {
 s/ +$//;
 my @t=split(/\t/);
 unless (@t>8) { print $_; next }
 my ($gid)=($t[8]=~m/gene_id "(MSTRG\.\d+)"/);
 if ($gid) {
   push(@prep, [$_, $gid]);
   my ($rn)=($t[8]=~m/ref_gene_id "([^"]+)/);
   #or for gene_name:
   #my ($rn)=($t[8]=~m/gene_name "([^"]+)/);
   if ($rn) {
     my $h=$g{$gid};
     if ($h) { $h->{$rn}=1 }
     else { $g{$gid}= { $rn=>1 } }
   }
 }
 else { print $_ }
}
my ($prevgid, $gadd);
foreach my $d (@prep) {
 my ($line, $gid)=@$d;
 if ($prevgid ne $gid) {
    $prevgid=$gid;
    $gadd='';
    if (my $gd=$g{$gid}) {
      $gadd='|'.join('|', (sort(keys(%$gd))));
    }
 }
 if ($gadd) {
   $line=~s/gene_id "MSTRG\.\d+/gene_id "$gid$gadd/;
   # 捕获原始transcript_id的后缀
   my ($transcript_suffix) = $line =~ /transcript_id "[^"]+\.(\w+)"/;
   # 只复制MSTRG ID部分到transcript_id，并保持原始后缀
   my ($base_gid) = $gid =~ /(MSTRG\.\d+)/;
   $line=~s/transcript_id "[^"]+"/transcript_id "$base_gid$gadd.$transcript_suffix"/;
 }
 print $line;
}

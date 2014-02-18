#!perl -w

my $ref = shift(@ARGV) or die "Error: no reference file name specified\n";

my $ref2 = shift(@ARGV) or die "Error: no new reference file name specified\n";


system(qq{rename_reference $ref $ref2 chrom_names.csv chrom_seqs.csv});

system(qq{bwa index $ref2});

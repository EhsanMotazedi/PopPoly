#!/usr/bin/env perl
# A simple perl wrapper for PopPoly to be able to specify the contig name and maximum genotype missing rate
# for haplotyping.
# Useful when multiple contigs are present in the bam and VCF files, which causes error for haplotyping.
# This wrapper extracts the desired contig from the bam and VCF files, extracts SNPs and breaks down complex varinats
# in the VCF file, filters out SNPs based on their genotype missing rate and runs PopPoly with the passed input options.
# This script should be stored in the same folder as PopPoly main.py and break_vcf.sh.
# Written by Ehsan Motazedi, Wageningen University & Research, 20-07-2018.
# Last modified: 15-08-2018 

use strict;
use warnings;
use File::Temp qw/ tempfile tempdir /;
use Cwd 'abs_path';
use Scalar::Util qw(looks_like_number);

sub signal_handler{
    die "\nCaught a signal: $!\n";
}

my @PATH = split /\//, abs_path($0);
my $PopPoly_Path = join '/', @PATH[0..$#PATH-1];
my $hlp_msg = "\nA simple perl wrapper for PopPoly, so that the contig name and the maximum genotype missing rate\nof each SNP can be specified for haplotyping. Useful when multiple contigs are present in the BAM\nand VCF files, which are not handled directly by PopPoly. Also, useful to exclude SNPs from phasings\nbased upon their genotype missing rate in the VCF file.\n\nThis wrapper extracts the desired contig from the bam and VCF files, performs SNP filtering and runs PopPoly\nwith the passed input options. To run this script it should be stored in the same folder as PopPoly main.py.\n\nInput parameters:\n\n\tparameters passed to PopPoly (see ./main.py -h)\n\n\t-h, --help\tprint this help message and exit.\n\t-c, --contig\tname of the contig to be phased.\n\t--maxMISS\tmaximum tolerable genotype missing rate for each SNP to be included in phasing.\n\nWritten by Ehsan Motazedi, Wageningen University & Research, 20-07-2018.";

$SIG{INT}  = \&signal_handler;
$SIG{TERM} = \&signal_handler;

my %options_set;
my @to_remove;
my @iftrue_opts = ("t","top","skip","filter","impute","redose");
my @toget_opts = ("c","contig","a", "aggressive", "e", "error", "r", "rho", "k", "kappa", "w", "warmup", "v", "verror", "P1", "P2", "mmq", "mbq", "maxIS", "maxMISS");
my $n=-1;

my ($fh1, $contig_bam) = tempfile(TEMPLATE => "tmp_bamXXXX", DIR=>"/tmp", SUFFIX => ".bam"); # tmp file to store contig reads if contig is specified
close($fh1);
my ($fh2, $contig_vcf) = tempfile(TEMPLATE => "tmp_vcfXXXX", DIR=>"/tmp", SUFFIX => ".vcf"); # tmp file to store contig variants if contig is specified
close($fh2);
my ($fh3, $contig_vcf_broken) = tempfile(TEMPLATE => "tmp_vcf_brokenXXXX", DIR=>"/tmp", SUFFIX => ".vcf"); # tmp file to store contig variants if contig is specified
close($fh3);

while ($n < $#ARGV){ #separate positional & optional arguments
	$n+=1;
	if (grep {/^-/} $ARGV[$n]){
		(my $option= $ARGV[$n]) =~ s/^-+//;
		if ($option eq "h" or $option eq "help"){
			print $hlp_msg."\n";
			exit 0;
		}
		if (grep {$_ eq $option} @iftrue_opts){
			$options_set{$option} = "";
			push @to_remove, $n;
		} elsif (grep {$_ eq $option} @toget_opts){
			if ($n == $#ARGV) {die "ERROR: No value passed to $ARGV[$n]!\n"}
			if (grep {/^-/} $ARGV[$n+1]) {die "ERROR: Invalid value $ARGV[$n+1] passed to $ARGV[$n]!\n"}
			$options_set{$option} = $ARGV[$n+1];
			push @to_remove, ($n,$n+1);
			$n+=1
		} else {die "ERROR: Invalid option $ARGV[$n]!\n"}
	}
}

for (reverse @to_remove){ # remove optional arguments (now saved in %options_set) from the input to get positional arguments
	splice @ARGV, $_, 1
}

if ($#ARGV == -1){ # check correct number of positional arguments
	die "ERROR: No bam file, VCF file and output name is given!\n"
} elsif ($#ARGV < 1){
	die "ERROR: No VCF file and output name is given!\n"
} elsif ($#ARGV < 2){
	die "ERROR: No output name is given!\n"
} elsif ($#ARGV > 2){
	die "ERROR: Too many positional arguments given!\n"
}

for (my $i=0; $i<14; $i+=2){ # check doubly given value getting options
	if (exists $options_set{$toget_opts[$i]} and exists $options_set{$toget_opts[$i+1]}) {die "ERROR: Only one of -$toget_opts[$i] and --$toget_opts[$i+1] allowed!\n"}
}
for (my $i=0; $i<2; $i+=2){ # check doubly given set true options 
	if (exists $options_set{$iftrue_opts[$i]} and exists $options_set{$iftrue_opts[$i+1]}) {die "ERROR: Only one of -$iftrue_opts[$i] and --$iftrue_opts[$i+1] allowed!\n"}
}

if (exists $options_set{maxMISS}){ # check the given SNP missing rate
	if (not looks_like_number($options_set{maxMISS})){
		die "ERROR: the SNP missing rate must be a number!\n"
	} elsif ($options_set{maxMISS}>1 or $options_set{maxMISS}<0) {
		die "ERROR: the SNP missing rate must be between 0 and 1 (0<=maxMISS<=1)!\n"
	} else {
		$options_set{maxMISS}+=1e-06
	}
} else {
	print STDERR "WARNING: No maximum missing rate is specified and hence no filtering of SNPs will be performed!\n";
	$options_set{maxMISS}=2
}

if (exists $options_set{contig} or exists $options_set{c}){ # check if contig name is given or not
	if (not exists $options_set{contig}){
		$options_set{contig}=$options_set{c};
		delete $options_set{c}
	}
}

(-f $ARGV[0]) or die "ERROR: could not find the bam file ${ARGV[0]}!\n";
(-f $ARGV[1]) or die "ERROR: could not find the VCF file ${ARGV[1]}!\n";

if (exists $options_set{contig}){  # extract contig reads from the original BAM fmile 
	if (! -f "$ARGV[0].bai"){
		print STDERR "WARNING: Cannot find the index of the input bam file: $ARGV[0].bai. samtools index is being run...\n";
		my $indx=system("samtools index $ARGV[0]");
		$indx ? die "ERROR: Indexing failed!\n$!\n" : print STDERR "Succeffully indexed the bam file...\n"
	}
	my $get_reads = system("samtools idxstats $ARGV[0]|grep -m 1 $options_set{contig} |awk '{print \"\\\"\"\$1\":0-\"\$2\"\\\"\"}'|xargs samtools view -bh $ARGV[0] > $contig_bam");
	!$get_reads or die "Cannot extract $options_set{contig} reads from $ARGV[0]: $!\n";
	$ARGV[0] = $contig_bam;
} else {
	print STDERR "WARNING: No contig name has been specified! The BAM and VCF files must contain only one contig!\n"
}

if ($options_set{maxMISS}<2 or exists $options_set{contig}) { # filter VCF file if contig and/or maxMISS are given
	my $ctg_found = 0; 
	open($fh2,">",$contig_vcf);  # extract contig variants from the original VCF file 
	open(my $fh_vcf, "<", $ARGV[1]) or die "Cannot open $ARGV[1]: $!\n";
	SEARCH: {
		while (my $line = <$fh_vcf>){ 
			chomp $line;
			if (grep {/^#/} $line){
				printf $fh2 "%s\n", $line;
			} else {
				my $vcf_ctg = (split /\t/, $line)[0];
				if (!(exists $options_set{contig}) or $vcf_ctg eq $options_set{contig}) {
					$ctg_found+=1;
					my $missing = 0;
					if ($options_set{maxMISS}<2) {  # filter SNPs using maxMISS rate
						my @genos = split /\t/, $line;
						my @non_missing_genos = grep {/([0-9](\/|\|))+/} map {(split /:/, $_)[0]} @genos[9..$#genos];
						$missing = ($#{genos}-9-$#{non_missing_genos})/($#{genos}+1-9);
					}
					if ($missing < $options_set{maxMISS}){
						printf $fh2 "%s\n", $line
					}
				} elsif ($ctg_found>0) { # avoid reading the rest of the VCF file after finding the desired contig
					last SEARCH;
				}
			}
		}
	}
	close($fh_vcf);
	close($fh2);
}

my $break_vcf = !system("${PopPoly_Path}/break_vcf.sh ${contig_vcf} ${contig_vcf_broken}"); # break complex variants and throw out indels using break_vcf.py
$break_vcf or die "Couldn't break the complex variants and throw out indels from the VCF file: $!\n";
$ARGV[1] = $contig_vcf_broken;

if (exists $options_set{contig}) {delete $options_set{contig}} # not passed to PopPoly main
if (exists $options_set{maxMISS}) {delete $options_set{maxMISS}} # not passed to PopPoly main
	
my %new_options_set;
for (keys %options_set){
	if (length($_)==1){
		$new_options_set{"-$_"}=$options_set{$_}
	} else {
		$new_options_set{"--$_"}=$options_set{$_}
	}
}
undef %options_set;
my $phasing = system((join " ","${PopPoly_Path}/main.py", @ARGV)." ".(join " ", map { $new_options_set{$_} eq "" ? "$_" : "$_ $new_options_set{$_}" } (keys %new_options_set))); 
!$phasing or die "ERROR: failed to run PopPoly: $!\n";

END {
    unlink "${contig_bam}" if (-f "${contig_bam}");
    unlink "${contig_vcf}" if (-f "${contig_vcf}");
    unlink "${contig_vcf_broken}" if (-f "${contig_vcf_broken}");
}

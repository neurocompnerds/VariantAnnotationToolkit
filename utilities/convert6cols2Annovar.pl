#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

our $REVISION = '$Revision: d0d72aa1e51ca9f3dc6b0611a4972956d6e93dec $';
our $DATE =	'$Date: 2014-07-14 03:13:01 -0700 (Mon, 14 Jul 2014) $';  
our $AUTHOR =	'$Author: Kai Wang <kai@openbioinformatics.org> $';

our ($verbose, $help, $man);
our ($variantfile);
our ($outfile, $format, $includeinfo, $snpqual, $snppvalue, $coverage, $maxcoverage, $chrmt, $altcov, $allelicfrac, $fraction, $species, 
	$filterword, $confraction, $allallele, $withzyg, $comment, $allsample, $genoqual, $varqual, $dbsnpfile, $withfreq, $withfilter, $seqdir, $inssize, $delsize, $subsize, $genefile, $splicing_threshold, $score);

our %iupac = (R=>'AG', Y=>'CT', S=>'CG', W=>'AT', K=>'GT', M=>'AC', A=>'AA', C=>'CC', G=>'GG', T=>'TT', B=>'CGT', D=>'AGT', H=>'ACT', V=>'ACG', N=>'ACGT', '.'=>'-', '-'=>'-'); ### <<< FOR 5500SOLiD LifeScope ( S=>'GC' is replaced by S=>'CG')
our %iupacrev = reverse %iupac; ### <<< FOR 5500SOLiD LifeScope

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'format=s'=>\$format, 'includeinfo'=>\$includeinfo,
	'snpqual=f'=>\$snpqual, 'snppvalue=f'=>\$snppvalue, 'coverage=i'=>\$coverage, 'maxcoverage=i'=>\$maxcoverage, 'chrmt=s'=>\$chrmt, 
	'fraction=f'=>\$fraction, 'altcov=i'=>\$altcov, 'allelicfrac'=>\$allelicfrac,
	'species'=>\$species, 'filter=s'=>\$filterword, 'confraction=f'=>\$confraction, 'allallele!'=>\$allallele, 'withzyg'=>\$withzyg,
	'comment'=>\$comment, 'allsample'=>\$allsample, 'genoqual=f'=>\$genoqual, 'varqual=f'=>\$varqual, 'dbsnpfile=s'=>\$dbsnpfile, 'withfreq'=>\$withfreq,
	'withfilter'=>\$withfilter, 'seqdir=s'=>\$seqdir, 'inssize=i'=>\$inssize, 'delsize=i'=>\$delsize, 'subsize=i'=>\$subsize, 'genefile=s'=>\$genefile,
	'splicing_threshold=i'=>\$splicing_threshold) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");

($variantfile) = @ARGV;

$chrmt ||= 'M';

if (defined $outfile) {
	open (STDOUT, ">$outfile") or die "Error: cannot write to output file $outfile: $!\n";
}

my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv) = qw/0 0 0 0 0 0 0 0 0/;

if ($variantfile eq 'stdin') {
	*VAR = *STDIN;
} elsif ($variantfile =~ m/^\.gz$/) {
	open (VAR, "gunzip -c $variantfile |") or die "Error: cannot read from STDIN uncompressing variant file $variantfile: $!\n";
} else {
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
}

while (<VAR>) {
	$countline++;
		
	if ($comment) {
		m/^#/ and print and next;
	} else {
		m/^#/ and next;		#skip comment lines
	}
	
	s/[\r\n]+$//;		#delete trailing new lines
	my $otherinfo = $_;	#this is the complete line (when -includeinfo is set, the entire line will be included in output file)
	
	my @field=split(/\t/);
	@field >=6 or die "Error: Likley not even a bastard VCF (at least 6 tab-delimited fields expected): <$_>\n";
	my ($chr, $start, $ID, $ref_allele, $mut_allele, $score) = @field;
	my ($end);
	my ($mut_allele2);
		
	#sometimes the alleles are not in the same case
	#chr1    1869771 1869774 actc    aCTctc          43.5    13      INDEL;DP=13;AF1=0.5;CI95=0.5,0.5;DP4=0,4,4,0;MQ=37;PV4=0.029,0.45,1,0.46
	$ref_allele = uc $ref_allele;
	$mut_allele = uc $mut_allele;
		
	#if ($ID eq '.' || $ID =~ /^rs/) {		#per MISHIMA, Hiroyuki suggestion (vcf4's third column (ID column) are not always ".")
	#	$end = $start;				#this block is commented out on 2011feb19
	#}
		
	if ($mut_allele eq '.') {			#no variant call was made at this position
		next;
	}
		
	if ($mut_allele =~ m/([^,]+),([\w,]+)/) {	#there could be more than two alternative alleles
		$mut_allele = $1;
		$mut_allele2 = $2;
	}
		
	if(length($ref_allele)==1 && length($mut_allele)==1) {  	### output snv
		if ($ref_allele =~ m/[^ACGTacgt]/ ) {
			print STDERR "WARNING: invalid allele record found in this file (ACGT expected): <$ref_allele> and <$mut_allele> in line <$_>\n";
			$ref_allele = 0;
			}
		if ( $mut_allele =~ m/[^ACGTacgt]/) {
			print STDERR "WARNING: invalid allele record found in this file (ACGT expected): <$ref_allele> and <$mut_allele> in line <$_>\n";
			$mut_allele = 0;
		}
			
		if ($includeinfo) {
			print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\t", $score, "\n";
		} else {
			print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\n";
		}
		
		if ($allallele) {
			if ($mut_allele2) {
				my @mut_allele2 = split (/,/, $mut_allele2);
				for my $i (0 .. @mut_allele2-1) {
					if ($includeinfo) {
						print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele2[$i], "\t", $score, "\n";
					} else {
						print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele2[$i], "\n";
					}
				}
			}
		}
		
		$countsnp++;
	} elsif (length($ref_allele) > 1 || length($mut_allele) > 1) {  ### output indel
			
		if(length($ref_allele) > length ($mut_allele)) { 		# deletion or block substitution
			my $head = substr($ref_allele, 0, length ($mut_allele));
			if ($head eq $mut_allele) {
				print $chr,"\t";
				print $start+length($head),"\t";
				print $start+length($ref_allele)-1,"\t";
				
				my $ref_allele1 = substr ($ref_allele, length ($mut_allele));
				print $ref_allele1,"\t";
				print "-";
			} else {
				print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele;
			}
		} elsif(length($mut_allele) >= length ($ref_allele)) { 		# insertion or block substitution
			my $head = substr ($mut_allele, 0, length ($ref_allele));
			if ($head eq $ref_allele) {
				print $chr,"\t";	
				print $start+length($ref_allele)-1,"\t";
				print $start+length($ref_allele)-1,"\t";
				
				$mut_allele = substr ($mut_allele, length ($ref_allele));
				print "-\t";
				print $mut_allele;
			} else {
				print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele;
			}
		}
			
		
		if ($includeinfo) {
			 print "\t", $score;
		} 
		print "\n";
		$countindel++;

		#do the same thing again, exactly like above, except that we work on second mutation;
		#in the future, consider rewrite this paragraph to make the code more elegant	
		if ($allallele and $mut_allele2) {
			my @mut_allele2 = split (/,/, $mut_allele2);
			for my $mut_allele2 (@mut_allele2) {
				if(length($ref_allele) > length ($mut_allele2)) { 		# deletion or block substitution
					my $head = substr($ref_allele, 0, length ($mut_allele2));
					if ($head eq $mut_allele2) {
						print $chr,"\t";
						print $start+length($head),"\t";
						print $start+length($ref_allele)-1,"\t";
						
						my $ref_allele1 = substr ($ref_allele, length ($mut_allele2));
						print $ref_allele1,"\t";
						print "-";
					} else {
						print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele2;
					}
				} elsif(length($mut_allele2) > length ($ref_allele)) { 		# insertion or block substitution
					my $head = substr ($mut_allele2, 0, length ($ref_allele));
					if ($head eq $ref_allele) {
						print $chr,"\t";	
						print $start+length($ref_allele)-1,"\t";
						print $start+length($ref_allele)-1,"\t";
						
						$mut_allele2 = substr ($mut_allele2, length ($ref_allele));
						print "-\t";
						print $mut_allele2;
					} else {
						print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele2;
					}
				} else {		#identical length of alleles
					print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele2;
					}
										
				if ($includeinfo) {
					print "\t", $score;
					} 
				print "\n";

			}
		}
	}
	$countvar++;
}


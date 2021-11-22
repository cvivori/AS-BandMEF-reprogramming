#!/bin/perl
$help="\n\tUsage: perl VTS_run_compare_all_couples.pl  path_INCL  min_dPSI  min_range  species DUPS/TRIPS/SINGLE  cond_1 ... cond_n\n
	\t1) path_INCL = path to inclusion table (output of vast-tools combine)
	\t2) min_dPSI = minimum delta PSI (see vast-tools compare)
	\t3) min_range = minimum distance between the ranges of groups (see vast-tools compare)
	\t4) species = species Hsa/Hs2/Mmu/Mm2/Dre.. (see vast-tools align)
	\t5) DUPS/TRIPS/MERGE = run vast-tools compare on duplicates, triplicates or single replicates 
	\t\t* choose single replicates to run on the merge-output!
	\t6+) conditions... = names of conditions 
	\t\t* excluded the suffixes _1 _2 _3 of replicates!

** Example run in the cluster: 
perl ~/Scripts/VTS_run_compare_all_couples.pl /users/jvalcarcel/cvivori/Transdifferentiation_BtoMac_GregJose/VASTTOOLS_v222/vast_out/INCLUSION_LEVELS_FULL-Hsa12-hg19.tab 10 5 Hsa DUPS BLaER_00h BLaER_24h BLaER_48h BLaER_96h\n\n died";


if ($#ARGV<6){
die $help; 
}

$INCL=shift; #$ARGV[0]
$dPSI=shift;
$range=shift;
$species=shift; 
$reps=shift;
die $help unless ($reps =~ /DUPS/ || $reps =~ /TRIPS/ || $reps =~ /SINGLE/) ;

$num = $#ARGV +1;
print "... Considering $num conditions ($reps):\n\t";
foreach $argnum (0 .. $#ARGV) {
    push @CONDs, $ARGV[$argnum];
}
print (join("\n\t",@CONDs));
print "\n";


if ($reps eq TRIPS) {

foreach $cond_1 (@CONDs) {
  foreach $cond_2 (@CONDs) {

	$c1r1= ($cond_1 . "_1");
	$c2r1= ($cond_2 . "_1");
	$c1r2= ($cond_1 . "_2");
	$c2r2= ($cond_2 . "_2");
	$c1r3= ($cond_1 . "_3");
	$c2r3= ($cond_2 . "_3");
	# next if $saw{$c1r3};
	# next if $saw{$c2r3};
	next if $saw{$cond_1}{$cond_2}>0;
	next if $saw{$cond_2}{$cond_1}>0;
	next if $cond_2 eq $cond_1;
	$saw{$cond_1}{$cond_2}++;
	$saw{$cond_2}{$cond_1}++;

print `qsub_job 6 6 6 -N vtCmp vast-tools compare $INCL -a $c1r1,$c1r2,$c1r3 -b $c2r1,$c2r2,$c2r3 --min_dPSI $dPSI --min_range $range --GO -sp $species --print_dPSI --print_sets`;

	}
	}

} elsif ($reps eq DUPS) {

foreach $cond_1 (@CONDs) {
  foreach $cond_2 (@CONDs) {

	$c1r1= ($cond_1 . "_1");
	$c2r1= ($cond_2 . "_1");
	$c1r2= ($cond_1 . "_2");
	$c2r2= ($cond_2 . "_2");
	next if $saw{$cond_1}{$cond_2}>0;
	next if $saw{$cond_2}{$cond_1}>0;
	next if $cond_2 eq $cond_1;
	$saw{$cond_1}{$cond_2}++;
	$saw{$cond_2}{$cond_1}++;
 
print `qsub_job 6 6 6 -N vtCmp vast-tools compare $INCL -a $c1r1,$c1r2 -b $c2r1,$c2r2 --min_dPSI $dPSI --min_range $range --GO -sp $species --print_dPSI --print_sets`;

	}
	}

} elsif ($reps eq SINGLE) {

foreach $cond_1 (@CONDs) {
  foreach $cond_2 (@CONDs) {

	$c1 = $cond_1;
	$c2= $cond_2;
	next if $saw{$cond_1}{$cond_2}>0;
	next if $saw{$cond_2}{$cond_1}>0;
	next if $cond_2 eq $cond_1;
	$saw{$cond_1}{$cond_2}++;
	$saw{$cond_2}{$cond_1}++;
 
print `qsub_job 6 6 6 -N vtCmp vast-tools compare $INCL -a $c1 -b $c2 --min_dPSI $dPSI --min_range $range --GO -sp $species --print_dPSI --print_sets`;

	}
	}

} else {die $help};


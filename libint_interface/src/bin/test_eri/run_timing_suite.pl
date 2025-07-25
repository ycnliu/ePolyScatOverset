#!/usr/bin/perl

my $lmax = 3;
my $ltotal_max = 10;
my $shellquartet_set = 1;

my @MakeVars = ("CXXGEN='icpc -m64 -std=c++0x' CXXGENFLAGS='-O3 -xCORE-AVX-I -ansi-alias'",
		"CXXGEN='clang++ -m64 -std=c++0x' CXXGENFLAGS='-O3 -march=native -mavx'");
my @UnrollingThresholds = ("4294967295", "0");
my @VectorLengths = ("1");
my @VectorizationMethods = ("0");
my @DoCSE = ("0", "1");

foreach my $makevars (@MakeVars) {

printf STDOUT "------------------------------------------------------------\n";
printf STDOUT "makevars: $makevars\n";
printf STDOUT "------------------------------------------------------------\n";

for($la=0; $la<=$lmax; ++$la) {
  for($lb=0; $lb<=$lmax; ++$lb) {
    for($lc=0; $lc<=$lmax; ++$lc) {
      for($ld=0; $ld<=$lmax; ++$ld) {
      
        my $ltotal = $la+$lb+$lc+$ld;
        next if $ltotal == 0 || $ltotal >= $ltotal_max;

        my $skip;
        if ($shellquartet_set == 1) {
          $skip = ! ShellQuartetPredicateStandard($la,$lb,$lc,$ld);
        }
        if ($shellquartet_set == 2) {
          $skip = ! ShellQuartetPredicateOrca($la,$lb,$lc,$ld);
        }
        
        next if $skip;

        foreach my $unroll_thresh (@UnrollingThresholds) {
          foreach my $veclength (@VectorLengths) {
            foreach my $vecmethod (@VectorizationMethods) {
              foreach my $cse (@DoCSE) {
                
                my $args = "$la $lb $lc $ld $unroll_thresh $veclength $vecmethod $cse";
                my $errcod = system("./test_eri.pl $args --makevars=\"$makevars\" 1>/dev/null 2>/dev/null");
                if ($errcod == 0) {
                  my $flop_output = `./time_eri 1 | grep nflops`;
                  $flop_output =~ s/\n//;
                  my $timing_output = `./run_time_eri.pl | grep FLOP`;
                  $timing_output =~ s/\n//;
                  printf STDOUT "$args $flop_output $timing_output\n";
	        }
	        else {
	          printf STDOUT "$args failed\n";
	        }
	      }
	    }
	  }
	}
	
      }
    }
  }
}

}

exit 0;

sub ShellQuartetPredicateStandard {
  my $la = shift;
  my $lb = shift;
  my $lc = shift;
  my $ld = shift;
  
  return $la >= $lb &&
         $lc >= $ld &&
         ($la+$lb <= $lc+$ld);
}

sub ShellQuartetPredicateOrca {
  my $la = shift;
  my $lb = shift;
  my $lc = shift;
  my $ld = shift;
  
  return $la <= $lb &&
         $lc <= $ld &&
         ($la < $lc || ($la == $lc && $lb <= $ld));
}

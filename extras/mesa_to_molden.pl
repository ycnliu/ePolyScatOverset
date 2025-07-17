#!/usr/bin/perl

$debug = 0;
$sentinel = 1;
while(<>){
    if ($sentinel){
        print $_ if $debug;
        if (/center\s+(\d+)\s+(\-*\d+\.\d+)\s+(\-*\d+\.\d+)\s+(\-*\d+\.\d+)\s+(\d+\.\d+)/){
            print "READING CENTER\n" if $debug;
            push(@centerxs, $2);
            push(@centerys, $3);
            push(@centerzs, $4);
            push(@charges, $5);
        }
        if (/basis fcn \s*(\d+).*are \s*(\d+)\s+(\d+)/){
            print "READING BASIS FCN FIRST/LAST\n" if $debug;
            push(@nfirsts, $2);
            push(@nlasts, $3);
        }
        if (/\s+(\d+)\s+(\d+)\s+(\d+)\s+(\-*\d+\.\d+)\s+(\-*\d+\.\d+)\s+(\-*\d+\.\d+)\s+(\-*\d+\.\d+)\s+(\-*\d+\.\d+)/){
            print "READING BASIS FCN INFO INC. EXP/COEF\n" if $debug;
            push(@ls, $1);
            push(@ms, $2);
            push(@ns, $3);
            push(@fcenterxs, $4);
            push(@fcenterys, $5);
            push(@fcenterzs, $6);
            push(@exps, $7);
            push(@coefs, $8);
            push(@dome, $1 == $1+$2+$3);
        }
    }
    if (/indexed/){
        print "ENDING FILE I/O\n" if $debug;
        $sentinel=0;
    }
}

@LSymbols = ("s", "p", "d", "f", "g", "h");
foreach $fcn(0..$#nfirsts){
    if (($fcenterxs[$nfirsts[$fcn]-1] != $fcenterxs[$nfirsts[$fcn]-2])||
        ($fcenterys[$nfirsts[$fcn]-1] != $fcenterys[$nfirsts[$fcn]-2])||
        ($fcenterzs[$nfirsts[$fcn]-1] != $fcenterzs[$nfirsts[$fcn]-2])){
            print "NEW CENTER FCN ", $fcn+1, "\n" if $debug;
            foreach $center(0..$#centerxs){
                if (($fcenterxs[$nfirsts[$fcn]-1] == $centerxs[$center])&&
                    ($fcenterys[$nfirsts[$fcn]-1] == $centerys[$center])&&
                    ($fcenterzs[$nfirsts[$fcn]-1] == $centerzs[$center])){
                    print "\tDETERMINED CENTER WAS CENTER ", $center+1, "\n" if $debug;
                    printf("\n %2i%2i\n",$center+1,0) if $dome[$nfirsts[$fcn]-1];
                    last;
                }
            }
        }
    foreach $prim($nfirsts[$fcn]-1..$nlasts[$fcn]-1){
        if ($prim == $nfirsts[$fcn]-1){
            print "FIRST PRIMITIVE\n" if $debug;
            $L = $ls[$prim] + $ms[$prim] + $ns[$prim];
            printf(" %-4s%2i%5.2f\n",$LSymbols[$L],$nlasts[$fcn]-$nfirsts[$fcn]+1,1.0) if $dome[$nfirsts[$fcn]-1];
        }
        printf("%18.10E%18.10E\n",$exps[$prim],$coefs[$prim]) if $dome[$nfirsts[$fcn]-1];
    }
}

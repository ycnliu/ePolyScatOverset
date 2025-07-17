#!/usr/bin/perl

while(<>){
    if (/\[MO\]/){
        $on=1;
    }
    if ($on){
        if (/^(\d+)/){
            if ($1==1){
                $count=$count+1;
            }
            $gcmd = "grep \"\\W$1\\W\" tmmop.out";
            $scmd = "sed -n '$count p'";
            # print "COMMANDS $gcmd GREP $scmd SED\n";
            $other = `$gcmd | $scmd`;
            # print "! $count $_ !! $other\n";
            print $other;
        }
        else{
            print;
        }
    }
}

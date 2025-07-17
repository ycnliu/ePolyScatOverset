#!/usr/bin/perl

$index = 0;
while(<>){
    if (/PRINT.*DBG/){
        # print "THERE HAS BEEN AN ERROR IN ONE OF THE RAY PRINTING STATEMENTS\n";
        # exit;
    }
    elsif (/PRINTING/){
        @split = split;
        $proc = $split[0];
        $dz = $split[5];
        if (not exists $dz_hash{$dz}){
            $idz = $idz + 1;
            $dz_hash{$dz} = $idz;
        }
        $rray = $split[8];
        # if ($rray == 1){
        #     print "\n !! SHOUT !!\n $_ \n\n";
        # }
        if (not $r_hash{$rray}){
            $ir = $ir + 1;
            $r_hash{$rray} = $ir;
        }
        # print "    DBG_DZH $dz_hash{$dz} $idz    ";
        if (($dz_hash{$dz} == 1) and ($idz > 1) and ($dz != $dz_last)){
        # if (($rray < $rray_last[$proc]) and (exists $dz_hash{$dz})){
        # print;
        # if ($rray < $rray_last){
        # print " SHOUT RRAY \n";
        # print " EXISTS? $dz ", exists $dz_hash{$dz}, "\n";
            $index = $index + 1;
        }
        $ugr = $split[11];
        $ugi = $split[12];
        $sgr = $split[15];
        $sgi = $split[16];
        # print;
        # print "$proc $dz_hash{$dz} $index $rray ($ugr, $ugi) ($sgr, $sgi)\n";
        $ugrmatrix[$dz_hash{$dz}-1][$proc][$index][$r_hash{$rray}-1] = $ugr;
        $ugimatrix[$dz_hash{$dz}-1][$proc][$index][$r_hash{$rray}-1] = $ugi;
        $sgrmatrix[$dz_hash{$dz}-1][$proc][$index][$r_hash{$rray}-1] = $sgr;
        $sgimatrix[$dz_hash{$dz}-1][$proc][$index][$r_hash{$rray}-1] = $sgi;
        $dz_last = $dz;
        # $rray_last = $rray;
    }
}

# $,="\n";
# print keys %r_hash;
# print "\n-----\n";
# print sort {$r_hash{$a} <=> $r_hash{$b}} keys %r_hash;
# print "\n-----\n";

foreach $iindex(0..$index){
    foreach $idz(0..$#ugrmatrix){
        print "# INDEX $iindex IDZ $idz \n";
        foreach $rval(sort {$r_hash{$a} <=> $r_hash{$b}} keys %r_hash){
            # print " !!! HI RVAL $rval !!! \n ";
            printf("%21.14E ", $rval);
            $ugrsum = 0.0;
            $ugisum = 0.0;
            $sgrsum = 0.0;
            $sgisum = 0.0;
            $ir = $r_hash{$rval} - 1;
            foreach $iproc(0..$#{$ugrmatrix[$idz]}){
                # printf("%.14E %.14E ",$ugmatrix[$idz][$iproc][$iindex][$ir],$sgmatrix[$idz][$iproc][$iindex][$ir]);
                $ugrsum = $ugrsum + $ugrmatrix[$idz][$iproc][$iindex][$ir];
                $ugisum = $ugisum + $ugimatrix[$idz][$iproc][$iindex][$ir];
                $sgrsum = $sgrsum + $sgrmatrix[$idz][$iproc][$iindex][$ir];
                $sgisum = $sgisum + $sgimatrix[$idz][$iproc][$iindex][$ir];
            }
            $totalr = $ugrsum + $sgrsum;
            $totali = $ugisum + $sgisum;
            if ($sgrsum**2+$sgisum**2 == 0.0){
                printf("%21.14E %21.14E %21s %21s %21.14E %21.14E\n",$ugrsum,$ugisum,"x.xE0","x.xE0",$totalr,$totali);
            }
            else{
                printf("%21.14E %21.14E %21.14E %21.14E %21.14E %21.14E\n",$ugrsum,$ugisum,$sgrsum,$sgisum,$totalr,$totali);
            }
        }
        print "\n\n";
    } 
    print "\n\n";
}




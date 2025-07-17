#!/usr/bin/perl

$DUMPFIL = @ARGV[0];
$PTSFIL = @ARGV[1];

open(PTSFIL,"<$PTSFIL");
while(<PTSFIL>){
    @split = split;
    $ptlist{$split[0]} = $split[1];
    # print "0 $split[0] 1 $split[1]\n";
}

foreach $pw(1..20){
    $col = 8 + ($pw-1)%5 + 1;
    $row = 5*int(($pw-1)/5) + 1;
    $grepcmd = "grep -n pw\\ *$row\\  $DUMPFIL | sed s/:/\\ /g | grep -v eg | sed s/[0-9]\\\\.[0-9]*\\\\-[0-9]*/0\\.00000000E+00/g";
    $awkcmd = "awk '{print \$6,\$1,\$$col}'";
    $sortcmd = "sort -n | awk '{print \$1, \$3}'";
    # print "$col/$row G $grepcmd A $awkcmd S $sortcmd\n";
    @rlist = ` $grepcmd | $awkcmd | $sortcmd `;
    $ptlast = -1;
    # print @rlist;
    printf("\n\n\n# PW $pw \{ $col/$row CMD $grepcmd | $awkcmd | $sortcmd \} \n");
#    foreach $r(@rlist){
#      print "R $r\n\n";
#      @split = split(/\s+/,$r);
#      $pt = $split[0];
#      $val = $split[1];
#      if ($pt != $ptlast){
#          print "!! NEWPOINT P $pt V $val!!\n";
#          print "!! NEWPOINT R $ptlist{$pt} V $val!!\n";
#      } else {
#          print "! OLDPOINT P $pt V $val !\n";
#      }
#      $ptlast = $pt;
#    }
    foreach $r(@rlist){
      @split = split(/\s+/,$r);
      $pt = $split[0];
      $val = $split[1];
      if ($pt != $ptlast){
          printf("\n%20.12E %20.12E ",$ptlist{$pt},$val);
      } else {
          printf("%20.12E ",$val);
      }
      $ptlast = $pt;
    }
}

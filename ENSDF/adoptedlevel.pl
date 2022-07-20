#!/usr/bin/perl

@element = (
  "  ",
  "H ", "HE", "LI", "BE", "B ", "C ", "N ", "O ", "F ", "NE",
  "NA", "MG", "AL", "SI", "P ", "S ", "CL", "AR", "K ", "CA",
  "SC", "TI", "V ", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
  "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y ", "ZR",
  "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN",
  "SB", "TE", "I ", "XE", "CS", "BA", "LA", "CE", "PR", "ND",
  "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB",
  "LU", "HF", "TA", "W ", "RE", "OS", "IR", "PT", "AU", "HG",
  "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH",
  "PA", "U ", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM",
  "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS",
  "RG", "CN", "NH", "FL", "MC", "LV", "TS", "OG");


open(FILE,@ARGV[0]);
while(<FILE>){
    if(/   ADOPTED LEVELS/){
        $id   = $_;
        $anum = substr($id,0,3);
        $symb = substr($id,3,2);
        for($i=1 ; $i<=$#element ; $i++){
            if($element[$i] eq $symb){
                $znum = $i;
            }
        }
#       print $znum," ",$anum," ",$symb,"\n";

        $k = 0;
        while(<FILE>){
            if(/^     /){ last; }
            $line[$k++] = $_;
        }
        $file = sprintf("ENSDF%03d%03d.dat",$znum,$anum);
        open(OUT,"> $file");
        print OUT $id;
        for($i=0 ; $i<$k ; $i++){
            print OUT $line[$i];
        }
        close(OUT);
    }
}

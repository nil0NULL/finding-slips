#!perl

##########################   Copyright and Licenses   ##########################
#                                                                              #
#   Copyright(C) 2021 Xudong He                                                #
#                                                                              #
#   This program is free software: you can redistribute it and/or modify       #
#   it under the terms of the GNU General Public License as published by       #
#   the Free Software Foundation, either version 3 of the License, or          #
#   at your option) any later version.                                         #
#                                                                              #
#   This program is distributed in the hope that it will be useful,            #
#   WITHOUT ANY WARRANTY; without even the implied warranty of                 #
#   or FITNESS FOR A PARTICULAR PURPOSE.  See the                              #
#   General Public License for more details.                                   #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.     #
#                                                                              #
################################################################################

################################################################################
#                                                                              #
#   Author: hxd_yi (Xudong He)                                                 #
#   E-mail: hxd_yi@163.com                                                     #
#   Description: This program is used to transform the cmf format to the       #
#   isparallel program input format.                                           #
#   Usage: RunMatScript.sh format [ cmf files ] | [ dirs ]                     #
#                                                                              #
################################################################################

use strict;
use warnings;
use English;
use Getopt::Long;
use MaterialsScript qw(:all);

use constant ABSZERO => 1e-6;

our $maxhkl = 10;
our $scalar_min = 0.8;
our $vdW_minus = 1.0;

our @vdwradii = (0.00, 
1.10,                                                                                                                                                                                     1.40, 
1.82, 1.53,                                                                                                                                                 1.92, 1.70, 1.55, 1.52, 1.47, 1.54, 
2.27, 1.73,                                                                                                                                                 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 
2.75, 2.31,                                                                                     2.15, 2.11, 2.07, 2.06, 2.05, 2.04, 2.00, 1.97, 1.96, 2.01, 1.87, 2.11, 1.85, 1.90, 1.85, 2.02, 
3.03, 2.49,                                                                                     2.32, 2.23, 2.18, 2.17, 2.16, 2.13, 2.10, 2.10, 2.11, 2.18, 1.93, 2.17, 2.06, 2.06, 1.98, 2.16, 
3.43, 2.68, 2.43, 2.42, 2.40, 2.39, 2.38, 2.36, 2.35, 2.34, 2.33, 2.31, 2.30, 2.29, 2.27, 2.26, 2.24, 2.23, 2.22, 2.18, 2.16, 2.16, 2.13, 2.13, 2.14, 2.23, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 
3.48, 2.83, 2.47, 2.45, 2.43, 2.39, 2.41, 2.43, 2.44, 2.45, 2.44, 2.45, 2.45, 2.45, 2.46, 2.46, 2.46);

our @covradii = (0.00, 
0.32,                                                                                                                                                                                     0.37, 
1.30, 0.99,                                                                                                                                                 0.84, 0.75, 0.71, 0.64, 0.60, 0.62,
1.60, 1.40,                                                                                                                                                 1.24, 1.14, 1.09, 1.04, 1.00, 1.01, 
2.00, 1.74,                                                                                     1.59, 1.48, 1.44, 1.30, 1.29, 1.24, 1.18, 1.17, 1.22, 1.20, 1.23, 1.20, 1.20, 1.18, 1.17, 1.16, 
2.15, 1.90,                                                                                     1.76, 1.64, 1.56, 1.46, 1.38, 1.36, 1.34, 1.30, 1.36, 1.40, 1.42, 1.40, 1.40, 1.37, 1.36, 1.36, 
2.38, 2.06, 1.94, 1.84, 1.90, 1.88, 1.86, 1.85, 1.83, 1.82, 1.81, 1.80, 1.79, 1.77, 1.77, 1.78, 1.74, 1.64, 1.58, 1.50, 1.41, 1.36, 1.32, 1.30, 1.30, 1.32, 1.44, 1.45, 1.50, 1.42, 1.48, 1.46, 
2.42, 2.11, 2.01, 1.90, 1.84, 1.80, 1.83, 1.80, 1.73, 1.68, 1.68, 1.68, 1.65, 1.67, 1.73, 1.76, 1.61);

if( ! scalar(@ARGV) ){
    $ARGV[0] = "."
}

for my $doc (@ARGV){
    print "$doc\n";
    if( -f $doc ){
        &xxx("$doc") if( -f $doc && $doc =~ /cmf$/i );
    }elsif( -d "$doc" ){
        opendir(my $dh, $doc) or die "Can't Open Direction $doc $! \n";
        while( readdir $dh ){
            #print "$doc/$_\n" if( -f $_ && /cmf$/i );
            if( $_ eq '.' || $_ eq '..' ){next;}
            &xxx("$doc/$_") if( -f "$doc/$_" && /cmf$/i );
        }
        closedir $dh;
    }
}

sub xxx
{
    my $inf = shift;
    my $prefix = $inf =~ s/\.cmf$//ir;
    #print "inf = $inf\n";
    #print "prefix = $prefix\n";
    my $doc = Documents->New("${prefix}.xsd") or die "Can't Create ${prefix}.xsd\n";
    $doc->CopyFrom($Documents{$inf}) or die "Can't Open $inf\n";
    print "Processing $inf\n";

    my $cell = $doc->UnitCell;
    my $output = Documents->New("${prefix}.txt") or die "Can't Create ${prefix}.txt\n";
    
    $output->Append(sprintf "%15.9f %15.9f %15.9f\n", $_->X, $_->Y, $_->Z) for ($cell->Lattice3D->VectorA, $cell->Lattice3D->VectorB, $cell->Lattice3D->VectorC);
    $output->Append(sprintf "%5d\n", $cell->Atoms->Count);

    eval {$doc->AssignMolecules;};
    if( $@ ) {
        $output->Append( &print_atoms($_) ) for @{$cell->Atoms};
    }else{
	$doc->GenerateLatticeDisplay(["PeriodicDisplayType"=>"Default"]);
        $cell = $doc->UnitCell;
        for my $mol_id ( 1 .. $cell->Molecules->Count ) {
            $output->Append( &print_atoms($_, $mol_id) ) for @{$cell->Molecules($mol_id - 1)->Atoms};
        }
    }
    #printf "%5d\n", $doc->UnitCell->Molecules->Count;

    #$output->Append(sprintf "%3d %-2s %9.6f %9.6f %9.6f %7.4f %5.2f %5.2f\n", $_->AtomicNumber, $_->ElementSymbol, $_->FractionalXYZ->X, $_->FractionalXYZ->Y, $_->FractionalXYZ->Z, $_->Mass, $covradii[$_->AtomicNumber], $vdwradii[$_->AtomicNumber]) for @{$cell->Atoms};
    $doc->Delete;
    return 0;
}

sub print_atoms
{
    my $x = shift;
    my $ans = sprintf "%3d %-2s %9.5f %9.5f %9.5f %8.4f %5.2f %5.2f", $x->AtomicNumber, $x->ElementSymbol, $x->FractionalXYZ->X, $x->FractionalXYZ->Y, $x->FractionalXYZ->Z, $x->Mass, $covradii[$x->AtomicNumber], $vdwradii[$x->AtomicNumber];

    my $mol_id;
    my $str;
    if( scalar @_ ) {
        $mol_id = shift;
        $str = sprintf " %5d\n", $mol_id;
    }else{
        $str = sprintf "\n";
    }
    return $ans . $str;
}


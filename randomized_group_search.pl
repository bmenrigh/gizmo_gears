#!/usr/bin/perl

use strict;
use warnings;

use Math::Trig;

require 'random_rc4-16.pl';


my $PI = 3.14159265358979323846264;

my $ORBIT_MAX = 32;


my %grips = ();
#add_grip('A', 3, 2.4, -1, 0);
#add_grip('B', 3, 2.4, 1, 0);

my $pcounter = 0;
my %points = ();


print 'Read("orbit_description.gap");', "\n";
for (my $i = 0; $i < 10000; $i++) {
    random_puzzle();
    sample_random_points();
}


sub random_puzzle {

    %grips = ();

    add_grip('A',
	     3 + rand_range_short_mask(2, 0x0F),
	     1.0 + rand_decimal() * 4.0,
	     -1, 0);

    add_grip('B',
	     3 + rand_range_short_mask(2, 0x0F),
	     1.0 + rand_decimal() * 4.0,
	     1, 0);

    print sprintf('Print("N1=%d; R1=%0.9f; N2=%d, R2=%.09f\n");',
		  $grips{'A'}{'N'}, $grips{'A'}{'R'},
		  $grips{'B'}{'N'}, $grips{'B'}{'R'}), "\n";

}


sub sample_random_points {

    for (my $i = 0; $i < 200; $i++) {
	find_orbit(random_point());
    }
}


sub find_orbit {
    my $x = shift;
    my $y = shift;

    %points = ();
    $pcounter = 0;

    add_point($x, $y, 0);

    my $added;
    my $d = 0;
    do {
	$added = do_depth($d);
	$d++;

	if ($added < 0) {
	    #warn 'Maxed out on points', "\n";
	    return;
	}

	#warn 'At depth ', $d, ' found ', $added, ' states', "\n";
    }
    while ($added > 0);

    return if ($pcounter < 6);

    # Find the cycles for each grip
    my @griplist;
    foreach my $grip (keys %grips) {
	my @permlist;
	foreach my $point (keys %points) {
	    my $newpoint =
		point_to_string(turn_grip_point($grip, $points{$point}{'X'},
						$points{$point}{'Y'}));

	    $permlist[$points{$point}{'C'} - 1] = $points{$newpoint}{'C'};
	}

	my $gripname = sprintf('grip_%s', $grip);
	print sprintf('%s := AsPermutation(Transformation([%s]));;',
		      $gripname, join(', ', @permlist)), "\n";

	push @griplist, $gripname;
    }

    #my $pname = point_to_string($x, $y);
    #$pname =~ tr/.-/pn/;
    my $pname = 'point';

    # Checking that each of these is a permutation shouldn't be needed
    # but the point rounding that is done can mess up and merge close points
    # into one which messes up the permutation
    foreach my $grip (@griplist) {
	print sprintf('if (IsPerm(%s) = true) then', $grip), "\n";
    }

    print sprintf('g_%s := Group([%s]);;',
		  $pname, join(', ', @griplist)), "\n";

    print sprintf('DescribeInterestingOrbits(g_%s);;', $pname), "\n";

    # Close each IsPerm() if statement
    foreach my $grip (@griplist) {
	print 'fi;', "\n";
    }

}


sub do_depth {
    my $d = shift;

    my $added = 0;

    foreach my $point (grep {$points{$_}{'D'} == $d} keys %points) {

	foreach my $grip (keys %grips) {

	    my ($nx, $ny) = turn_grip_point($grip, $points{$point}{'X'},
					    $points{$point}{'Y'});

	    my $pcount = add_point($nx, $ny, $d + 1);

	    if ($pcount > 0) {
		$added++;

		if ($pcount > $ORBIT_MAX) {
		    return -1;
		}
	    }
	}
    }

    return $added;
}



sub add_grip {
    my $name = shift;
    my $n = shift;
    my $r = shift;
    my $x = shift;
    my $y = shift;

    my $theta = (2.0 * $PI) / ($n * 1.0);

    if (exists $grips{$name}) {
	die 'Got duplicate grip name ', $name, "\n";
    }

    $grips{$name} = {()};
    $grips{$name}{'N'} = $n * 1.0;
    $grips{$name}{'R'} = $r * 1.0;
    $grips{$name}{'X'} = $x * 1.0;
    $grips{$name}{'Y'} = $y * 1.0;
    $grips{$name}{'theta'} = $theta;
    $grips{$name}{'costheta'} = cos($theta);
    $grips{$name}{'sintheta'} = sin($theta);
}


sub dist {
    my ($x1, $y1, $x2, $y2) = @_;

    return sqrt((($x2 - $x1) ** 2) + (($y2 - $y1) ** 2));
}


sub ingrip {
    my $grip = shift;
    my $x = shift;
    my $y = shift;

    if (dist($x, $y,
	     $grips{$grip}{'X'}, $grips{$grip}{'Y'}) < $grips{$grip}{'R'}) {
	return 1;
    }
    else {
	return 0;
    }
}


sub random_point {

    my @griplist = keys %grips;

    fisher_yates(\@griplist);

    return random_point_in_grip($griplist[0]);
}


sub random_point_in_grip {
    my $grip = shift;

    my $r = rand_decimal() * $grips{$grip}{'R'};
    my $t = rand_decimal() * (2.0 * $PI);

    my $x = $grips{$grip}{'X'} + $r * cos($t);
    my $y = $grips{$grip}{'Y'} + $r * sin($t);

    return ($x, $y);
}


sub turn_grip_point {
    my $grip = shift;
    my $x = shift;
    my $y = shift;

    # If the point is in the grip
    if (ingrip($grip, $x, $y) == 1) {
	# Shift the point over to be centered on the origin
	$x = $x - $grips{$grip}{'X'};
	$y = $y - $grips{$grip}{'Y'};

	# Rotation in 2D about the origin via rec coords:
	# x' = x*cos(t) - y*sin(t)
	# y' = x*sin(t) + y*cos(t)
	my $nx = $x * $grips{$grip}{'costheta'} -
	    $y * $grips{$grip}{'sintheta'};
	my $ny = $x * $grips{$grip}{'sintheta'} +
	    $y * $grips{$grip}{'costheta'};

	# Shift the point back
	$x = $nx + $grips{$grip}{'X'};
	$y = $ny + $grips{$grip}{'Y'};
    }

    return ($x, $y);
}


sub add_point {
    my $x = shift;
    my $y = shift;
    my $d = shift; # found depth

    my $point_str = point_to_string($x, $y);

    if (exists $points{$point_str}) {
	return 0; # Point already added
    }

    $pcounter++;

    $points{$point_str} = {()};
    $points{$point_str}{'X'} = $x;
    $points{$point_str}{'Y'} = $y;
    $points{$point_str}{'D'} = $d;
    $points{$point_str}{'C'} = $pcounter;

    return $pcounter;
}


sub point_to_string {
    my $x = shift;
    my $y = shift;

    return sprintf('%.09f_%.09f', $x, $y);
}

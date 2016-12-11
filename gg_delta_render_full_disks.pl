#!/usr/bin/perl

use strict;
use warnings;

use GD;
use Math::Trig;
use Math::Round;
use File::Glob;
use File::Copy;

my $parallelism = 4;
my $children = 0;

my $PI = 3.14159265358979323846264;

my $param_n = 5;
#my $param_r = sqrt((7 + sqrt(5)) / 2);
my $param_r = sqrt((7 + sqrt(5)) / 2);
#my $param_r_txt = 'sqrt((7 + sqrt(5)) / 2)';
my $param_r_txt = 'sqrt((7 + sqrt(5)) / 2)';

my $aturn = 4;
my $bturn = 1;

# Negative because turns for this puzzle are clockwise
# but positive rotations are counter-clockwise
my $atheta = (($aturn * -1.0) / ($param_n * 1.0)) * (2.0 * $PI);
my $btheta = (($bturn * -1.0) / ($param_n * 1.0)) * (2.0 * $PI);
my $ntheta = ((1.0) / ($param_n * 1.0)) * (2.0 * $PI);
my $cosatheta = cos($atheta);
my $sinatheta = sin($atheta);
my $cosbtheta = cos($btheta);
my $sinbtheta = sin($btheta);


# === These settings work for showing the disks
#my $goalw = 512;
#my $scalef = ($goalw * 1.0) / ((2.0 * $param_r) + 2.0);

#my $w = int(((2.0 * $param_r) + 2.0) * $scalef);
#my $h = int(2.0 * $param_r * $scalef);

#my ($xmin, $xmax) = (-1.0 - $param_r, 1.0 + $param_r);
#my ($ymin, $ymax) = (-1.0 * $param_r, 1.0 * $param_r);
#===

# === These settings work for showing the wedge
my $wheight = sqrt(($param_r ** 2) - 1);
my $wwidth = $param_r - 1;

my ($xmin, $xmax) = (-1.0 * $wwidth, 1.0 * $wwidth);
my ($ymin, $ymax) = (-1.0 * $wheight, 1.0 * $wheight);

my $h = 512;
my $w = int(($wwidth / $wheight) * $h);
# ===

# === These settings work for "zooming" into a smaller region
#my $w = 2048;
#my $h = 2048;
#
#my ($xmin, $xmax) = (-0.2, 0.2);
#my ($ymin, $ymax) = (-0.2, 0.2);
#===


my $ih = $h; # Actual image height accounting for legend

my $delta_color = 0; # Are colors based on delta, not absolute order?
my $smooth_delta = 1; # Should color be smooth across ratio of a/b near 1?
my $border_color = 1; # Are colors based on border distance, not order?
my $blend_border = 1; # Should border pixels get blended with black?
my $add_color_legend = 1;
my $legend_pad = 16;
my $legend_height = 32;
if ($add_color_legend == 1) {
    $ih = $h + $legend_pad + $legend_height;
}

# We have to substract 1 from the width and hight here because
# the point 0, 0 is actually in the center of a pixel
# which means we have each pixel stick out .5 to the top, bottom, left, and
# right of the disks / wedge.  This means we have 1 extra pixel worth of
# height and width we must account for
my $pwidth = (($xmax - $xmin) / (($w - 1.0) * 1.0));
my $pheight = (($ymax - $ymin) / (($h - 1.0) * 1.0));
my $pradius = sqrt(($pwidth ** 2.0) + ($pheight ** 2.0));

my $border_innerr = $param_r - (2.0 * $pradius);
my $border_outerr = $param_r + (2.0 * $pradius);

# Sampling controls
my $aa_samp = 16;
#my $aa_ord_cutoff = 200 * 1 * 1000;
my $border_samples = 1024;
my $samp_disk_points = 0;
my $samp_only_border = 0;
#my $missing_data_cutoff = 2 * 1000 * 1000; # Not used for delta A/B code
my $max_point_sample_cutoff = 16 * 1024 * 1024;
my $start_point_sample_cutoff = 64 * 1024;
my $min_wq_batch_size = 4 * $parallelism;
my $start_wq_batch_size = 4096 * $parallelism;
my $sample_round_factor = 2; # How much to change these limits each pass

my $point_sample_cutoff = $start_point_sample_cutoff;
my $wq_batch_size = $start_wq_batch_size;


my $lastsave = time();
my $saverate = 10 * 60; # 10 minutes

GD::Image->trueColor(1); # turn on true color
my $cimg = GD::Image->new(1, 1, 1); # preallocate some colors

my @white = (255, 255, 255);
my $white_idx = $cimg->colorAllocate(@white);
my @black = (0, 0, 0);
my $black_idx = $cimg->colorAllocate(@black);
my $cmin_idx;
my $cmax_idx;
if ($delta_color == 0) {
    $cmin_idx = $cimg->colorAllocate(val_to_rgb(0.0));
    $cmax_idx = $cimg->colorAllocate(val_to_rgb(0.0));
}
else {
    $cmin_idx = $cimg->colorAllocate(val_to_rgb(-1.0));
    $cmax_idx = $cimg->colorAllocate(val_to_rgb(1.0));
}

my $aa_log = 1; # Use the geometric mean rather than average

#my $usezbox = 1;
my $usezbox = 0;
my $sym180 = 1;
my $symmoves = 0; # Rotates points all the way around the disks
#my $sym180 = 0;

my $stats_pcount = 0;
my $stats_scount = 0;
my $stats_usable_scount = 0;
my $stats_sym180_scount = 0;
my $stats_symmoves_scount = 0;
my $stats_bogus_samples = 0;
my $stats_data_missing = 0;
my $stats_aa_pixels_needed = 0;
my $stats_aa_samples_needed = 0;
my $stats_border_samples_needed = 0;
my $stats_passes = 0;

my @igrid_order;
my @igrid_scount;

my @bpoints_froma_x;
my @bpoints_froma_y;
my @apoints_fromb_x;
my @apoints_fromb_y;
my $cent_min_r = 2.0 - $param_r;
my $cent_max_r = dist((cos($ntheta / 2.0) * $param_r) - 1.0,
		      sin($ntheta / 2.0) * $param_r,
		      1, 0);

$bpoints_froma_x[0] = 1.0;
$bpoints_froma_y[0] = 0.0;
$apoints_fromb_x[0] = -1.0;
$apoints_fromb_y[0] = 0.0;

for (my $i = 1; $i < $param_n; $i++) {

	my $x = $bpoints_froma_x[$i - 1] + 1.0;
	my $y = $bpoints_froma_y[$i - 1];

	my $nx = $x * cos($ntheta) - $y * sin($ntheta);
	my $ny = $x * sin($ntheta) + $y * cos($ntheta);

	$x = $nx - 1.0;
	$y = $ny;

	$bpoints_froma_x[$i] = $x;
	$bpoints_froma_y[$i] = $y;
	$apoints_fromb_x[$i] = $x * -1.0;
	$apoints_fromb_y[$i] = $y;
}


my $omin = 2000000000; # Large number. min should always be less than this ;-)
my $omax = 0;


my $OG_NAME = sprintf('gg_delta_state_ogrid_d%d_b%d_n%d_r%.8f_a%db%d_%dx%d_' .
		      '[%f,%f]_[%f,%f].ggs',
		      $delta_color, $border_color,
		      $param_n, $param_r, $aturn, $bturn, $w, $h,
		      $xmin, $xmax, $ymin, $ymax);

my $SC_NAME = sprintf('gg_delta_state_scount_d%d_b%d_n%d_r%.8f_a%db%d_%dx%d_' .
		      '[%f,%f]_[%f,%f].ggs',
		      $delta_color, $border_color,
		      $param_n, $param_r, $aturn, $bturn, $w, $h,
		      $xmin, $xmax, $ymin, $ymax);

my $IMG_NAME = sprintf('gg_delta_img_fulld_d%d_b%d_n%d_r%.8f_a%db%d_%dx%d_' .
		       '[%f,%f]_[%f,%f].png',
		       $delta_color, $border_color,
		       $param_n, $param_r,
		       $aturn, $bturn, $w, $h,
		       $xmin, $xmax, $ymin, $ymax);


# Used for point selection jitter
srand(time());

#####
read_state_files();
read_existing_points_files();
render_image();
#####

sub dist {
    my ($x1, $y1, $x2, $y2) = @_;

    return sqrt((($x2 - $x1) ** 2) + (($y2 - $y1) ** 2));
}


sub inwedge {
    my $x = shift;
    my $y = shift;

    if ((dist($x, $y, -1, 0) < $param_r) &&
	(dist($x, $y, 1, 0) < $param_r)) {
	return 1;
    }
    else {
	return 0;
    }
}


sub indisks {
    my $x = shift;
    my $y = shift;

    if ((dist($x, $y, -1, 0) < $param_r) ||
	(dist($x, $y, 1, 0) < $param_r)) {
	return 1;
    }
    else {
	return 0;
    }
}

sub incenter {
    my $x = shift;
    my $y = shift;

    my $adist = dist($x, $y, -1.0, 0);
    my $bdist = dist($x, $y, 1.0, 0);

    # Check in A
    if ($adist < $param_r) {
	return 1 if ($adist < $cent_min_r);
	return 0 if ($adist > $cent_max_r);

	for (my $i = 0; $i < $param_n; $i++) {
	    return 0 if (dist($x, $y,
			      $bpoints_froma_x[$i],
			      $bpoints_froma_y[$i]) < $param_r);
	}
	return 1;
    }
    elsif ($bdist < $param_r) { # Check B
	return 1 if ($bdist < $cent_min_r);
	return 0 if ($bdist > $cent_max_r);

	for (my $i = 0; $i < $param_n; $i++) {
	    return 0 if (dist($x, $y,
			      $apoints_fromb_x[$i],
			      $apoints_fromb_y[$i]) < $param_r);
	}
	return 1;
    }
    else {
	return 0;
    }
}


sub inborder {
    my $x = shift;
    my $y = shift;

    # A border pixel is in a 2 pixel radius of the outer rim of a disk
    my $dista = dist($x, $y, -1, 0);
    my $distb = dist($x, $y, 1, 0);

    if ((($dista >= $border_innerr) && ($dista <= $border_outerr)) ||
	(($distb >= $border_innerr) && ($distb <= $border_outerr))) {
	return 1;
    }
    else {
	return 0;
    }
}


sub val_to_rgb {
    my $v = shift;

    # $v must be in the range [0, 1]

    if (($v < -1.0) || ($v > 1.0)) {
	warn 'Bogus color value: ', $v, "\n";
	return (0, 0, 0);
    }

    my ($r, $g, $b);

    if (($delta_color == 0) || ($smooth_delta == 0)) {
	if ($v >= 0) {
	    $v *= ($PI / 2.0);

	    # Try to keep at least one color high valued and out of phase
	    # with the others so that colors stay saturated
	    $r = round(sin($v) * 255.0);
	    $g = round((1 - sin($v * 2.0)) * 255.0);
	    $b = round(cos($v) * 255.0);
	}
	else {
	    $v *= ($PI / -2.0);

	    #$r = int(sin($v * $PI) * 255.0);
	    #$g = int(((sin($v * ($PI * 2.0)) + 1.0) / 2) * 255.0);
	    #$b = int(((cos($v * ($PI * 2.0)) + 1.0) / 2) * 255.0);

	    $r = round((1 - sin($v)) * 255.0);
	    $g = round(sin($v * 2.0) * 255.0);
	    $b = round((1 - cos($v)) * 255.0);
	}
    }
    else {
	# This color function needs to be smooth about 1 (val near 0)
	# To do this, the range [-1, 1] is mapped into [0, Pi/2] by splitting
	# up [0, 1] into [Pi / 4, 0] (backwards)
	# and [-1, 0] into [Pi / 4, Pi / 2]
	if ($v >= 0) {
	    $v = 1.0 - $v;
	    $v *= ($PI / 4.0);
	}
	else {
	    $v *= ($PI / -4.0);
	    $v += ($PI / 4.0);
	}

	# The wavefactor is how many "bumps" to add to each channel
	# The waveamp is how strong the wavefactor applies
	my $wavefactor = 5;
	my $waveamp = 0.1;

	# Scale the wavefactor to fit in Pi / 2
	$wavefactor = ($wavefactor * 4.0) + 1.0;

	$r = round((sin($v) * (1.0 - $waveamp) +
		    sin($v * $wavefactor) * $waveamp) * 255.0);
	$g = round(((1 - sin($v * 2.0)) * (1.0 - $waveamp) +
		    (1 - sin($v * 2.0 * $wavefactor)) * $waveamp) * 255.0);
	$b = round((cos($v) * (1.0 - $waveamp) +
		    cos($v * $wavefactor) * $waveamp) * 255.0);
    }

    return ($r, $g, $b);
}


sub order_to_val {
    my $o = shift;
    my $aaomin = shift;
    my $aaomax = shift;

    my $v;

    if ($delta_color == 0) {
	#$v = atan2($o - $aaomin, 1.0) / atan2($aaomax - $aaomin, 1.0);

	my $offset = 5;

	$v = (log(($o - $aaomin) + $offset) /
	      log(($aaomax - $aaomin) + $offset))
    }
    else {
	if ($o > 1) {
	    my $op = $o;
	    my $omaxp = $aaomax;
	    #$v = (log(log($op) + 1) / log(log($omaxp) + 1));
	    $v = atan2($op - 1.0, 1.0) / atan2($omaxp - 1.0, 1.0);
	}
	else {
	    my $on = 1.0 / $o;
	    my $omaxn = 1.0 / $aaomin;
	    #$v = (log(log($on) + 1) / log(log($omaxn) + 1));
	    $v = atan2($on - 1.0, 1.0) / atan2($omaxn - 1.0, 1.0);
	    $v *= -1.0;
	}
    }

    # Due to rounding errors this can be a tiny bit larger than 1 which is
    # a no-no
    $v = 1.0 if ($v > 1.0);
    $v = -1.0 if ($v < -1.0);

    return $v;
}


sub point_to_ixiy {
    my $x = shift;
    my $y = shift;

    if (($x + ($pwidth / 2.0) < $xmin) || ($x - ($pwidth / 2.0) > $xmax) ||
	($y + ($pheight / 2.0) < $ymin) || ($y - ($pheight / 2.0) > $ymax)) {
	return (-1, -1);
    }

    my $ix = int((($x - $xmin) / $pwidth) + 0.5);
    my $iy = int((($y - $ymin) / $pheight) + 0.5);

    return ($ix, $iy);
}


sub ixiy_to_point {
    my $ix = shift;
    my $iy = shift;
    my $ixoff = shift;
    my $iyoff = shift;

    if (($ix < 0) || ($ix >= $w) ||
	($iy < 0) || ($iy >= $h)) {
	return (-1, -1);
    }

    my $x = (($pwidth * ($ix + $ixoff)) + $xmin) - ($pwidth / 2.0);
    my $y = (($pheight * ($iy + $iyoff)) + $ymin) - ($pheight / 2.0);

    return ($x, $y);
}


sub is_border_pixel {
    my $ix = shift;
    my $iy = shift;

    # A border pixel has at least one corner in the disks
    # and at least one corner out of the disks
    my $incount = 0;
    my ($x, $y);

    # Check corner 1
    ($x, $y) = ixiy_to_point($ix, $iy, 0, 0);
    $incount += 1 if (indisks($x, $y) == 1);

    # Check corner 2
    ($x, $y) = ixiy_to_point($ix, $iy, 0, 1);
    $incount += 1 if (indisks($x, $y) == 1);

    # Check corner 3
    ($x, $y) = ixiy_to_point($ix, $iy, 1, 0);
    $incount += 1 if (indisks($x, $y) == 1);

    # Check corner 4
    ($x, $y) = ixiy_to_point($ix, $iy, 1, 1);
    $incount += 1 if (indisks($x, $y) == 1);

    if (($incount > 0) && ($incount < 4)) {
	return 1;
    }

    return 0;
}


sub border_blend_amount {
    my $ix = shift;
    my $iy = shift;

    # A border pixel is a square and not all of its area is in the disks
    # Some of the area is outside of the puzzle and should be black
    # so this function samples a bunch of points in the pixel
    # to see how much of the pixel's area is inside of the disks

    my $samplecount = 1024;

    my ($x, $y);

    my $incount = 0;
    for (my $i = 0; $i < $samplecount; $i++) {
	($x, $y) = ixiy_to_point($ix, $iy, rand(), rand());
	$incount += 1 if (indisks($x, $y) == 1);
    }

    return (($incount * 1.0) / ($samplecount * 1.0));
}


sub move_point {
    my $x = shift;
    my $y = shift;

    # Rotation in 2D about the origin via rec coords:
    # x' = x*cos(t) - y*sin(t)
    # y' = x*sin(t) + y*cos(t)

    # If the point is in circle A
    if (dist($x, $y, -1, 0) < $param_r) {
	# Shift the point over to be centered on the origin
	$x = $x + 1.0;

	my $nx = $x * $cosatheta - $y * $sinatheta;
	my $ny = $x * $sinatheta + $y * $cosatheta;

	$x = $nx - 1.0;
	$y = $ny;
    }

    # If the point is in circle B
    if (dist($x, $y, 1, 0) < $param_r) {
	# Shift the point over to be centered on the origin
	$x = $x - 1.0;

	my $nx = $x * $cosbtheta - $y * $sinbtheta;
	my $ny = $x * $sinbtheta + $y * $cosbtheta;

	$x = $nx + 1.0;
	$y = $ny;
    }

    return ($x, $y);
}


sub aa_order_point {
    my $ix = shift;
    my $iy = shift;

    my $aa_order;
    if ($aa_log == 1) {
	$aa_order = exp((1.0 * $igrid_order[$ix][$iy]) /
			(1.0 * $igrid_scount[$ix][$iy]));
    }
    else {
	$aa_order = ((1.0 * $igrid_order[$ix][$iy]) /
		     (1.0 * $igrid_scount[$ix][$iy]));
    }

    return $aa_order;
}


sub neigh_avg_order {

    my $ix = shift;
    my $iy = shift;

    my $ncount = 0;
    my $nosum = 0;

    for (my $cix = $ix - 1; $cix <= $ix + 1; $cix++) {
	next if (($cix < 0) || ($cix >= $w));

	for (my $ciy = $iy - 1; $ciy <= $iy + 1; $ciy++) {
	    next if (($ciy < 0) || ($ciy >= $h));

	    next if (($cix == $ix) && ($ciy == $iy));

	    next if ($igrid_scount[$cix][$ciy] == 0);

	    $nosum += aa_order_point($cix, $ciy);
	    $ncount++;
	}
    }

    if ($ncount > 0) {
	return (($nosum * 1.0) / ($ncount * 1.0));
    }
    else {
	return -1;
    }
}


sub read_points_file {
    my $fname = shift;

    warn 'Reading points from ', $fname, "\n";

    open(my $inpoints, '<', $fname) or die 'Unable to open ', $fname, ' for ',
    'reading: ', $!, ' ', $?, "\n";

    while (<$inpoints>) {
	chomp;

	my $line = $_;

	if ($line =~ m/^(-?\d+\.\d{4,16})\s(-?\d+\.\d{4,16})\s
                       (\d+)\s(\d+)\s(\d+\.\d{4,16})$/x) {
	    my ($px, $py, $ordera, $orderb, $toborder) = ($1, $2, $3, $4, $5);

	    next if ($ordera == 0);
	    next if ($orderb == 0);

	    my ($ordr, $mordr);
	    if ($border_color == 0) {
		if ($delta_color == 0) {
		    $ordr = ($ordera * 1.0) + ($orderb * 1.0);
		    $mordr = $ordr;
		}
		else {
		    $ordr = ($ordera * 1.0) / ($orderb * 1.0);
		    $mordr = ($orderb * 1.0) / ($ordera * 1.0);
		}
	    }
	    else {
		if ($delta_color == 1) {
		    # TODO: make this actually work
		    # It needs to somehow scale numbers very close
		    # to zero to be very close to 1
		    # and numbers not so close to 0 need to be close
		    # to zero or huge depending on if aturn > bturn

		    # Scale proximity to large number
		    my $recip = 1.0 + (1.0 / $toborder);
		    if ($ordera > $orderb) {
			$ordr = $recip;
			$mordr = (1.0 / $ordr);
		    }
		    else {
			$ordr = 1.0 / $recip;
			$mordr = $recip;
		    }
		}
		else {
		    $ordr = (1.0 / $toborder);
		    $mordr = $ordr;
		}
	    }

	    add_point($px, $py, $ordr);

	    if ($symmoves == 1) {
		my $nx = $px;
		my $ny = $py;

		for (my $i = 0; $i < (2 * $param_n); $i++) {
		    ($nx, $ny) = move_point($nx, $ny);

		    if (inwedge($nx, $ny) == 0) {
			$stats_symmoves_scount++;
			add_point($nx, $ny, $ordr);
		    }
		    else {
			last;
		    }
		}
	    }


	    if ($sym180 == 1) {
		$stats_sym180_scount++;
		add_point($px * -1.0, $py * -1.0, $mordr);

		if ($symmoves == 1) {
		    my $nx = $px * -1.0;
		    my $ny = $py * -1.0;

		    for (my $i = 0; $i < (2 * $param_n); $i++) {
			($nx, $ny) = move_point($nx, $ny);

			if (inwedge($nx, $ny) == 0) {
			    $stats_symmoves_scount++;
			    add_point($nx, $ny, $mordr);
			}
			else {
			    last;
			}
		    }
		}

	    }

	}
	else {
	    $stats_bogus_samples++;
	}
    }

    close($inpoints);

    move($fname, 'processed/' . $fname) or
	die 'Unable to move ', $fname, ' ', $?, ' : ', $!, "\n";
}


sub add_point {
    my $px = shift;
    my $py = shift;
    my $order = shift;

    $stats_scount++;

    my ($ix, $iy) = point_to_ixiy($px, $py);

    #warn 'Got point: ', $px, ' ', $py, "\n";

    return if (($ix == -1) || ($iy == -1));

    if (($ix < 0) || ($ix >= $w) ||
	($iy < 0) || ($iy >= $h)) {
	$stats_bogus_samples++;
	#warn 'Got bogus pixel coords: ', $ix, ' ', $iy, "\n";
	return;
    }

    if ($order > $omax) {
	$omax = $order;
    }
    if ($order < $omin) {
	$omin = $order;
    }

    $stats_usable_scount++;

    if ($aa_log == 1) {
	$igrid_order[$ix][$iy] += log($order * 1.0);
    }
    else {
	$igrid_order[$ix][$iy] += $order * 1.0;
    }
    $igrid_scount[$ix][$iy] += 1;

    #warn 'Sample count for ', $ix, ' ', $iy, ' : ',
    #$igrid_scount[$ix][$iy], "\n";
}


sub find_order_min_max {

    # Short cuircut function to make sure colors are normal
    # across multiple rendering of differnt parts of the image
    #return (1.0 / 11.0, 11.0);

    my $aaomin = 2000000000; # Large number. min should always be less ;-)
    my $aaomax = 0;

    for (my $ix = 0; $ix < $w; $ix++) {
	for (my $iy = 0; $iy < $h; $iy++) {

	    if ($igrid_scount[$ix][$iy] > 0) {

		my $aao = aa_order_point($ix, $iy);

		if ($aao > $aaomax) {
		    $aaomax = $aao;
		}
		if ($aao < $aaomin) {
		    $aaomin = $aao;
		}
	    }
	}
    }

    if ($delta_color == 1) {
	# Due to divide-by-zero issues with atan2 these need to be
	# less than or greater than 1.0
	if ($aaomin >= 1) {
	    $aaomin = 0.999;
	}
	if ($aaomax <= 1) {
	    $aaomax = 1.001;
	}
    }

    return ($aaomin, $aaomax);
}


sub get_point_sample_cmd {
    my $x = shift;
    my $y = shift;

    if ($usezbox == 0) {
	return sprintf('pointcycle_ab_delta(%d, %s, [%.08f, %.08f, 1]~, ' .
		       '%d, %d, %d)',
		       $param_n, $param_r_txt,
		       $x, $y, $aturn, $bturn,
		       $point_sample_cutoff);
    }
    else {
	#return sprintf('pointcycleorderabfastzbox' .
	#	       '([%.08f, %.08f, 1]~, %d, %d, %d)',
	#	       $x, $y, $aturn, $bturn,
	#	       $point_sample_cutoff);
	die 'Zbox currently unimplemented!', "\n";
    }
}


#open(my $outsamp, '>', $ARGV[2]) or die 'Unable to open ', $ARGV[2], ' for ',
#    'writing: ', $!, ' ', $?, "\n";


#sub write_point_sample {
#    my $x = shift;
#    my $y = shift;
#
#    print $outsamp get_point_sample_cmd($x, $y), "\n";
#}


sub do_image_pass {

    my @work_queue = ();
    my ($msamp, $asamp) = (0, 0);
    my ($aaomin, $aaomax) = find_order_min_max();

    $stats_passes++;

    $stats_data_missing = 0;
    $stats_aa_pixels_needed = 0;
    $stats_aa_samples_needed = 0;
    $stats_border_samples_needed = 0;

    GD::Image->trueColor(1);

    # 1 indicates true color
    my $img = GD::Image->new($w, $ih, 1);

    warn 'Starting pass through pixel grid...', "\n";

    for (my $ix = 0; $ix < $w; $ix++) {
	for (my $iy = 0; $iy < $h; $iy++) {

	    if (scalar @work_queue > $wq_batch_size) {
		do_work(\@work_queue);
		@work_queue = ();
	    }

	    my $scount = $igrid_scount[$ix][$iy];

	    my ($cx, $cy) = ixiy_to_point($ix, $iy, 0.5, 0.5);

	    my $pinwedge = inwedge($cx, $cy);
	    my $pindisks = indisks($cx, $cy);
	    my $pincenter = incenter($cx, $cy);
	    my $pinborder = inborder($cx, $cy);

	    my $isborderpixel = 0;
	    if (($pinborder == 1) && ($blend_border == 1)) {
		$isborderpixel = is_border_pixel($ix, $iy);
	    }

	    # Unsampled data
	    if ($scount == 0) {

		# If the pixel is in the wedge
		if (($pinwedge == 1) ||
		    (($samp_disk_points == 1) &&
		     ($pindisks == 1) &&
		     ($pincenter == 0))) {

		    $stats_data_missing++;

		    $msamp++;

		    $img->setPixel($ix, ($ih - 1) - $iy, $white_idx);

		    if ($pinborder == 1) {
			#print $outsamp '/* Missing border data */', "\n";
			for (my $i = 0; $i < $border_samples; $i++) {
			    $stats_border_samples_needed++;


			    push @work_queue, get_point_sample_cmd(
				ixiy_to_point($ix, $iy,
					      rand(1), rand(1)));

			}
		    }
		    elsif ($samp_only_border == 0) {

			# we could only sample pixels that have all
			# their neigbors filled in or whoes surronding
			# pixels have loweish order but neither of
			# those make much sense for the deta A/B code.

			#my $neigh_avg = neigh_avg_order($ix, $iy);
			#if ($neigh_avg != -1) {
			#if ($neigh_avg < $missing_data_cutoff) {
			    #print $outsamp '/* Missing data sample */', "\n";

			push @work_queue, get_point_sample_cmd(
			    ixiy_to_point($ix, $iy,
					  rand(1), rand(1)));

			#}
			#}
		    }
		}
		elsif (($samp_disk_points == 1) && ($pincenter == 1)) {
		    # A gets max
		    if (dist($cx, $cy, -1.0, 0.0) < $param_r) {
			$img->setPixel($ix, ($ih - 1) - $iy, $cmax_idx);
		    }
		    else {
			$img->setPixel($ix, ($ih - 1) - $iy, $cmin_idx);
		    }
		}
		else {
		    # Not in the wedge (or disks or we aren't sampling
		    # non-wedge points

		    $img->setPixel($ix, ($ih - 1) - $iy, $black_idx);
		}
	    } # end sample count is zero
	    else { # samle count is not zero
		my $aao = aa_order_point($ix, $iy);

		# Due to additional sample gathering the order may
		# actually increase
		if ($aao > $aaomax) {
		    $aaomax = $aao;
		}
		# Same for the min
		if ($aao < $aaomin) {
		    $aaomin = $aao;
		}

		$stats_pcount++;

		# See if we may want to re-sample this point for
		# better AA
		if (($pinwedge == 1) &&
		    (($pinborder == 1) || ($samp_only_border == 0))){
		    #if ($aao < $aa_ord_cutoff) { # doesnt make sense for delta
		    if ($scount < $aa_samp) {
			my $neigh_avg = neigh_avg_order($ix, $iy);

			if ($neigh_avg != -1) {

			    # It is possible that points collected
			    # could have reduced the order of
			    # neighboring pixels but we haven't
			    # visited one of those pixels yet so we
			    # still think the minimum order is too
			    # large so we must check that the neighbor
			    # average is greater than the min order

			    if (($neigh_avg < $aaomin) ||
				(abs(order_to_val($neigh_avg, $aaomin,
						  $aaomax) -
				     order_to_val($aao, $aaomin,
						  $aaomax)) >
				 (1.0 / 255.0))) {

				$stats_aa_pixels_needed++;

				$asamp++;

				#print $outsamp '/* AA re-sample */', "\n";
				for (my $i = $scount; $i < $aa_samp; $i++) {
				    $stats_aa_samples_needed++;

				    push @work_queue, get_point_sample_cmd(
					ixiy_to_point($ix, $iy,
						      rand(1), rand(1)));

				}
			    }
			}
		    }
		}
		#warn 'Got aa order: ', $aa_order, "\n";

		my @c = val_to_rgb(order_to_val($aao, $aaomin, $aaomax));

		# If this is a border pixel we should blend with black
		if (($blend_border == 1) &&
		    ($pinborder == 1) &&
		    ($isborderpixel == 1)) {

		    my $blend_v = border_blend_amount($ix, $iy);

		    $c[0] = round($c[0] * $blend_v);
		    $c[1] = round($c[1] * $blend_v);
		    $c[2] = round($c[2] * $blend_v);
		}

		my $cidx = $img->colorAllocate(@c);

		$img->setPixel($ix, ($ih - 1) - $iy, $cidx);
	    }

	} # End for IY
    } # End for IX

    if ($add_color_legend == 1) {
	for (my $ix = 0; $ix < $w; $ix++) {
	    for (my $iy = $h; $iy < $h + $legend_pad; $iy++) {
		$img->setPixel($ix, ($ih - 1) - $iy, $black_idx);
	    }
	}
	if ($delta_color == 0) {
	    for (my $ix = 0; $ix < $w; $ix++) {
		for (my $iy = $h + $legend_pad; $iy < $ih; $iy++) {
		    my $legend_idx = $cimg->
			colorAllocate(val_to_rgb((1.0 / ($w - 1)) * $ix));
		    $img->setPixel($ix, ($ih - 1) - $iy, $legend_idx);
		}
	    }
	}
	else {
	    for (my $ix = 0; $ix < $w; $ix++) {
		for (my $iy = $h + $legend_pad; $iy < $ih; $iy++) {
		    my $legend_idx = $cimg->
			colorAllocate(val_to_rgb(((-2.0 / ($w - 1)) *
						  $ix) + 1.0));
		    $img->setPixel($ix, ($ih - 1) - $iy, $legend_idx);
		}
	    }
	}
    }



    if (scalar @work_queue > 0) {
	do_work(\@work_queue);
	@work_queue = ();
    }

    #


    open(my $outimg, '>', $IMG_NAME) or
	die 'Unable to open ', $IMG_NAME, ' for writing: ', $!, ' ', $?, "\n";

    print $outimg $img->png();

    close($outimg);

    write_state_files();

    warn '==== Pass ', $stats_passes, ' stats: ====', "\n";
    warn 'Samples count: ', $stats_scount, "\n";
    warn 'Samples from 180 symmetry: ', $stats_sym180_scount, "\n";
    warn 'Samples from moves symmetry: ', $stats_symmoves_scount, "\n";
    warn 'Usable samples: ', $stats_usable_scount, "\n";
    warn 'Bogus samples: ', $stats_bogus_samples, "\n";
    warn 'Pixels with samples: ', $stats_pcount, "\n";
    warn 'Pixels without samples: ', $stats_data_missing, "\n";
    warn 'Pixels needing anti-aliasing: ', $stats_aa_pixels_needed, "\n";
    warn 'Anti-aliasing samples needed: ', $stats_aa_samples_needed, "\n";
    warn 'Border samples needed: ', $stats_border_samples_needed, "\n";
    warn 'Minimum anti-aliased order seen: ', $aaomin, "\n";
    warn 'Maximum anti-aliased order seen: ', $aaomax, "\n";
    warn 'Minimum order seen: ', $omin, "\n";
    warn 'Maximum order seen: ', $omax, "\n";

    return ($msamp, $asamp);
}


sub do_work {
    my $wq_ref = shift;

    try_save();

    warn 'Starting work queue batch...', "\n";

    my %pnames;

    my $wq_len = scalar @{$wq_ref};

    for (my $n = 0; $n < $parallelism; $n++) {

	my $sesid = get_session_id();

	my $W_NAME = sprintf('gg_delta_samples_n%d_r%.8f_a%db%d_%s.gp',
			     $param_n, $param_r, $aturn, $bturn, $sesid);

	my $P_NAME = sprintf('gg_delta_points_n%d_r%.8f_a%db%d_%s.txt',
			     $param_n, $param_r, $aturn, $bturn, $sesid);

	my $P_TEMP = sprintf('gg_delta_points_n%d_r%.8f_a%db%d_%s.temp',
			     $param_n, $param_r, $aturn, $bturn, $sesid);

	$pnames{$P_NAME} = 1;

	open(my $outsamp, '>', $W_NAME) or
	    die 'Unable to open ', $W_NAME, ' for writing: ', $!, ' ', $?, "\n";
	for (my $i = $n; $i < $wq_len; $i += $parallelism) {
	    print $outsamp $wq_ref->[$i], "\n";
	}
	print $outsamp '\\q', "\n";
	close($outsamp);

	my $cmd = 'gp -q gg_gen_points_full_param.gp ' . $W_NAME .
	    ' > ' . $P_TEMP . ' 2>/dev/null';

	if (fork() == 0) {

	    my $ret = `$cmd`;

	    unlink($W_NAME);
	    move($P_TEMP, $P_NAME);

	    exit(0);
	}
	else {
	    $children++;
	}
    }

    while ($children > 0) {
	wait();
	$children--;

	foreach my $P_NAME (keys %pnames) {
	    if (-e $P_NAME) {
		read_points_file($P_NAME);
		delete $pnames{$P_NAME};
	    }
	}
    }

    if (scalar keys %pnames > 0) {
	warn 'Some points files seem to have been left over!', "\n";
    }
}


sub read_existing_points_files {

    my $P_GLOB = sprintf('gg_delta_points_n%d_r%.8f_a%db%d_*.txt',
			 $param_n, $param_r, $aturn, $bturn);

    foreach my $P_NAME (sort glob $P_GLOB) {
	read_points_file($P_NAME);
    }
}


sub write_state_files {

    warn 'Starting to write state...', "\n";

    open(my $OGFH, '>', $OG_NAME) or
	die 'Unable to open ', $OG_NAME, ' ', $?, ' : ', $!, "\n";

    open(my $SCFH, '>', $SC_NAME) or
	die 'Unable to open ', $SC_NAME, ' ', $?, ' : ', $!, "\n";


    for (my $iy = 0; $iy < $h; $iy++) {
	my @ogrow;
	my @scrow;
	for (my $ix = 0; $ix < $w; $ix++) {
	    push @ogrow, sprintf('%.8f', $igrid_order[$ix][$iy]);
	    push @scrow, $igrid_scount[$ix][$iy];
	}
	print $OGFH join(' ', @ogrow), "\n";
	print $SCFH join(' ', @scrow), "\n";
    }

    close($OGFH);
    close($SCFH);

    warn 'State written...', "\n";
}


sub read_state_files {

    if ((-e $OG_NAME) &&
	(-e $SC_NAME)) {

	warn 'Reading state from files...', "\n";

	open(my $OGFH, '<', $OG_NAME) or
	    die 'Unable to open ', $OG_NAME, ' ', $?, ' : ', $!, "\n";

	open(my $SCFH, '<', $SC_NAME) or
	    die 'Unable to open ', $SC_NAME, ' ', $?, ' : ', $!, "\n";

	@igrid_order = ();
	@igrid_scount = ();
	for (my $i = 0; $i < $w; $i++) {
	    my @orcol = ((0) x $h);
	    my @sccol = ((0) x $h);

	    push @igrid_order, [@orcol];
	    push @igrid_scount, [@sccol];
	}

	#my $scount = $igrid_scount[$ix][$iy];

	my $iy = 0;
	while (<$OGFH>) {
	    chomp;

	    my $line = $_;

	    next if ($line =~ m/^\s*#/);

	    my @row = split(/\s+/, $line);

	    for (my $ix = 0; $ix < $w; $ix++) {
		$igrid_order[$ix][$iy] = $row[$ix] * 1.0;
	    }

	    $iy++;
	}

	$iy = 0;
	while (<$SCFH>) {
	    chomp;

	    my $line = $_;

	    next if ($line =~ m/^\s*#/);

	    my @row = split(/\s+/, $line);

	    for (my $ix = 0; $ix < $w; $ix++) {
		$igrid_scount[$ix][$iy] = int($row[$ix]);
	    }

	    $iy++;
	}

    }
    else {

	warn 'Could not find state files, starting with empty state...', "\n";

	@igrid_order = ();
	@igrid_scount = ();
	for (my $i = 0; $i < $w; $i++) {
	    my @orcol = ((0) x $h);
	    my @sccol = ((0) x $h);

	    push @igrid_order, [@orcol];
	    push @igrid_scount, [@sccol];
	}
    }
}


sub try_save {

    my $curtime = time();

    if ($curtime - $lastsave >= $saverate) {
	write_state_files();
	$lastsave = $curtime;
    }
}


sub get_session_id {

    open(URANDOM, '<', '/dev/urandom') or
	die 'Unable to open urandom: ', $?, "\n";

    my $rand;
    read(URANDOM, $rand, 16);

    close URANDOM;

    return unpack('H*', $rand);
}


sub render_image {

    my ($prev_msamp, $prev_asamp) = do_image_pass();

    my $done = 0;
    if (($prev_msamp == 0) && ($prev_asamp == 0)) {
	warn 'No work left to be done!', "\n";
	$done = 1;
    }

    my $no_aa_progress = 0;
    while ($done == 0) {

	# Adjust paramaters each round
	$point_sample_cutoff = round($point_sample_cutoff *
				     $sample_round_factor);

	if ($point_sample_cutoff > $max_point_sample_cutoff) {
	    $point_sample_cutoff = $max_point_sample_cutoff;
	}

	$wq_batch_size = round($wq_batch_size /
			       $sample_round_factor);

	if ($wq_batch_size < $min_wq_batch_size) {
	    $wq_batch_size = $min_wq_batch_size;
	}

	# Now do a pass with the new parameters
	my ($msamp, $asamp) = do_image_pass();

	if (($msamp == 0) && ($asamp == 0)) {
	    warn 'No work left to be done!', "\n";
	    $done = 1;
	    last;
	}

	if (($msamp > 0) && ($msamp >= $prev_msamp)) {
	    warn 'No sampling progress made!', "\n";
	    $done = 1;
	    last;
	}

	if (($msamp == 0) && ($asamp > ($prev_asamp / 2.0))) {
	    #warn 'No anti-aliasing progress made!', "\n";
	    $no_aa_progress++;
	}

	if ($no_aa_progress > 4) {
	    warn 'No anti-aliasing progress being made!', "\n";
	    $done = 1;
	    last;
	}

	($prev_msamp, $prev_asamp) = ($msamp, $asamp);
    }

}

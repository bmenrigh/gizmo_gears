#!/usr/bin/perl

use strict;
use warnings;

no warnings "portable";

my $rc4_16_state;
my $rc4_16_i;
my $rc4_16_j;


sub rc4_16_init {
    my $key = shift;
    my $keylen = length $key;

    $rc4_16_state = pack('n' x 65536, (0 .. 65535));
    $rc4_16_i = 0;
    $rc4_16_j = 0;

    for ($rc4_16_i = 0; $rc4_16_i <= 65535; $rc4_16_i++) {
        $rc4_16_j = ($rc4_16_j +
                     unpack('n', substr($rc4_16_state, $rc4_16_i * 2, 2)) +
                     unpack('n', substr($key, ($rc4_16_i * 2) % $keylen, 1) .
                            substr($key, (($rc4_16_i * 2) + 1) % $keylen)))
            & 0xFFFF;

        # Swap
        (substr($rc4_16_state, $rc4_16_i * 2, 2),
         substr($rc4_16_state, $rc4_16_j * 2, 2)) =
             (substr($rc4_16_state, $rc4_16_j * 2, 2),
              substr($rc4_16_state, $rc4_16_i * 2, 2));

    }

    $rc4_16_i = 0;
    $rc4_16_j = 0;


    # throw away 4x the state to reduce initial bias
    for(my $i = 0; $i < (4 * 65536); $i++) {
        my $junk = rc4_16_word();
    }
}

sub rc4_16_word {
    $rc4_16_i = ($rc4_16_i + 1) & 0xFFFF;
    $rc4_16_j = ($rc4_16_j +
                 unpack('n', substr($rc4_16_state, $rc4_16_i * 2, 2)))
        & 0xFFFF;

    # Swap
    (substr($rc4_16_state, $rc4_16_i * 2, 2),
     substr($rc4_16_state, $rc4_16_j * 2, 2)) =
         (substr($rc4_16_state, $rc4_16_j * 2, 2),
          substr($rc4_16_state, $rc4_16_i * 2, 2));

    return substr($rc4_16_state,
                  ((unpack('n', substr($rc4_16_state, $rc4_16_i * 2, 2)) +
                    unpack('n', substr($rc4_16_state, $rc4_16_j * 2, 2))) &
                   0xFFFF), 2);
}


sub rc4_16_scalar {

    return unpack('n', rc4_16_word());
}


sub random_qword {

    return (rc4_16_word() . rc4_16_word() . rc4_16_word() . rc4_16_word());

}


sub rand_range_mask {
    my $max = shift;
    my $mask = shift;

    my $rlonglong;
    my $r;

    do {
        $rlonglong = unpack('Q', random_qword());
        $r = $rlonglong & $mask;
    } while ($r > $max);


    return $r;
}


sub rand_range_short_mask {
    my $max = shift;
    my $mask = shift;

    my $r16;
    my $r;

    while ($mask >> 1 > $max) {
        $mask >>= 1;
    }

    do {
        $r16 = rc4_16_scalar();
        $r = $r16 & $mask;
    } while ($r > $max);

    return $r;
}


sub rand_range {
    my $max = shift;

    if ($max == 0) {
        return 0;
    }

    my $mask = 0xFFFFFFFFFFFFFFFF;
    while ($mask >> 1 > $max) {
        $mask >>= 1;
    }
    #print 'mask: ', $mask, "\n";

    my $rlonglong;
    my $r;

    do {
        $rlonglong = unpack('Q', random_qword());
        #print 'rlonglong: ', $rlonglong, "\n";
        $r = ($rlonglong & $mask);
        #print 'r: ', $r, "\n";
    } while ($r > $max);


    return $r;
}


sub rand_decimal {
    return ((1.0 * unpack('Q', random_qword())) / (2.0 ** 64.0));
}


sub fisher_yates {
    my $listref = shift;

    my $len = scalar @{$listref};

    for (my $i = ($len - 1); $i >= 1; $i--) {
        my $j = rand_range($i);
        ($listref->[$i], $listref->[$j]) = ($listref->[$j], $listref->[$i]);
    }
}


sub fisher_yates_short_mask {
    my $listref = shift;
    my $mask = shift;

    my $len = scalar @{$listref};

    for (my $i = ($len - 1); $i >= 1; $i--) {
        my $j = rand_range_short_mask($i, $mask);
        ($listref->[$i], $listref->[$j]) = ($listref->[$j], $listref->[$i]);
    }
}


#print 'Qword: 0x', unpack('H*', random_qword()), "\n";

rc4_16_init($$ . 'whatever' . time() . $< . $( . $0);


sub rc4_16_reinit {
    rc4_16_init($$ . 'whatever' . time() . $< . $( . $0 . join('_', @_));
}


#for(my $i = 0; $i < 100; $i++) {
#    print rand_range(100), "\n";
#}



1;

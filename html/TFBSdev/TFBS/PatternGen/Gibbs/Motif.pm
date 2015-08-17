package TFBS::PatternGen::Gibbs::Motif;
use vars qw(@ISA);
use strict;

use TFBS::Matrix::PFM;
use TFBS::PatternGen::Motif;
@ISA = qw(TFBS::PatternGen::Motif);


sub _calculate_PFM  {
    my $self = shift;
    my @PFM;
    foreach my $rowref ( @{$self->{'matrix'}} )  {
	my @PFMrow;
	foreach my $element (@$rowref) {
	    push @PFMrow, int($self->{'nr_hits'}*$element/100 + 0.5);
	}
	push @PFM, [@PFMrow];
    }
    return [@PFM];
}


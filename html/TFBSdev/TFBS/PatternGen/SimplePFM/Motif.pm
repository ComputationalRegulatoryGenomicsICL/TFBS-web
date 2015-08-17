package TFBS::PatternGen::Motif;
use vars qw(@ISA);
use strict;

use TFBS::Matrix::PFM;

@ISA = qw(Bio::Root::RootI);

sub new  {
    my ($caller, %args) = @_;
    my $self = bless {}, ref($caller) || $caller;
    $self->{'matrix'} = 
	$args{'-matrix'} 
              || $self->throw("No -matrix provided");
    $self->{'length'} = $args{'-length'} || scalar @{$self->{'matrix'}->[0]};
    $self->{'nr_hits'} = $args{'-nr_hits'} 
               || $self->throw("No -nr_hits provided.");
    $self->{'bg_probabilities'} = 
	$args{'-bg_probabilities'} || [0.25,0.25,0.25,0.25];
    return $self;
}

sub PFM  {
    my ($self, %args) = @_;
    return TFBS::Matrix::PFM->new (-name => "unknown",
				   -ID   => "unknown",
				   -class=> "unknown",
				   %args,
				   -matrix => $self->_calculate_PFM()
				   );
}

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


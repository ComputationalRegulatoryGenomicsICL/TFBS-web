package TFBS::_Iterator;

use vars '@ISA';
use strict;
use Carp;
@ISA = qw(Bio::Root::RootI);

#############################################################
# PUBLIC METHODS
#############################################################

sub new  {
    my ($caller, $arrayref, $sort_by, $reverse) = @_;
    my $class = ref $caller || $caller;
    my $self;
    if ($arrayref)  {
	$self = bless { _orig_array_ref     => [ @$arrayref ],
			_iterator_array_ref => [ @$arrayref ],
			_sort_by            => ($sort_by || undef),
			_reverse            => ($reverse || 0)
			},
		$class;
    }
    else  {
	croak("No no valid array ref for Iterator of ".
	      (ref($class)  || $class)." provided:");
    }
    
    $self->_sort()    if $sort_by;
    $self->_reverse() if $reverse;

    return $self;
}
				       


sub current {

}

sub reset  {
    my ($self) = @_;
    @{$self->{_iterator_array_ref}} = @{$self->{_orig_array_ref}};
    $self->_sort()    if $self->{'_sort_by'};
    $self->_reverse() if $self->{'reverse'};
    return $self;
}

sub next {
    my $self = shift;
    return shift @{$self->{_iterator_array_ref}};
}
#################################################################
# PRIVATE METHODS
#################################################################

sub _sort  {
    my ($self, $sort_by) = @_;
    $sort_by or $sort_by = $self->{_sort_by} or  $sort_by = 'name';

    # we can sort by name, start, end, score
    my %sort_fn = 
	(start => sub  { $a->start() <=> $b->start() 
			     ||
			  $a->pattern->name() cmp $b->pattern->{name}
		       },
	 end   => sub  { $a->end()   <=> $b->end()   
			     ||
			  $a->pattern->name() cmp $b->pattern->{name}
		       },
	 name  => sub  { $a->end()   <=> $b->end()   
			     ||
			  $a->pattern->name() cmp $b->pattern->{name}
		       },
	 score => sub {  $b->score()   <=> $a->score()
			     ||
			  $a->pattern->name() cmp $b->pattern->name()
		      }
	);
			 
    if (defined (my $sort_function = $sort_fn{lc $sort_by})) {
	$self->{'_iterator_array_ref'} =
	    [ sort $sort_function @{$self->{'_orig_array_ref'}} ];
    }
    else  {
	$self->throw("Cannot sort ".ref($self)." object by '$sort_by'.");
    }
}

sub _reverse {
    my $self = shift;
    $self->{'_iterator_array_ref'} = 
	[ reverse @{ $self->{'_iterator_array_ref'} } ];
}





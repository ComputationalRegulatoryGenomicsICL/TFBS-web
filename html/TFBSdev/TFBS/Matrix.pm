# TFBS module for TFBS::Matrix
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

TFBS::Matrix - base class for matrix patterns, containing methods common to all

=head1 DESCRIPTION

TFBS::Matrix is a base class consisting of universal constructor called by its subclasses (TFBS::Matrix::*), and matrix manipulation methods that are independent of the matrix type. It is not meant to be instantiated itself.

=head1 FEEDBACK

Please send bug reports and other comments to the author.

=head1 AUTHOR - Boris Lenhard

Boris Lenhard E<lt>Boris.Lenhard@cgb.ki.seE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are preceded with an underscore.

=cut



# The code begins HERE:

package TFBS::Matrix;
use vars '@ISA';

use PDL; # this dependancy has to be elimiated in the future versions
use Bio::Root::RootI;

use strict;

@ISA = qw(Bio::Root::RootI);

sub new  {
    my $class = shift;
    my %args = @_;
    my $self = bless {}, ref($class) || $class;

    # first figure out how it was called
    # we need (-dbh and (-ID or -name) for fetching it from a database
    #         or -matrix for direct matrix input
    
    unless (defined $args{-matrixtype}) {
	$self->throw("No matrix type defined. \n".
		     "You shoud probably not create ".
		     "TFBS::Matrix object directly"); 
    }
    if (my $dbobj= $args{-db}) {

	# This is an obsolete way of fetching a matrix from a database.
	# It will be discontinued in future releases
       
	if (defined $args{'-name'})  {
	    $self = 
		$dbobj->get_Matrix_by_name($args{-name},$args{-matrixtype}) 
	      or return undef;
	}
	elsif (defined $args{-ID}) { 
	    $self = $dbobj->get_Matrix_by_ID ($args{-ID}, $args{-matrixtype})
	      or return undef;

	}
	else  {
	    $self->throw("No matrix ID or name provided ".
				"alongside the db object.");
	}
    }
    elsif (defined $args{'-matrix'}) {
	$self->set_matrix($args{'-matrix'});
    }
    elsif (defined $args{'-matrixstring'}) {
	$self->set_matrix($args{'-matrixstring'});
    }
    elsif (defined $args{-matrixfile}) {
	$self->set_matrix(`cat $args{-matrixfile}`); # FIXME: TEMPORARY
    }
    else  {
	$self->throw("No matrix or db object provided.");
    }

    # Set the object data.
    # Parameters specified in constructor call override those
    # fetched from the database.

    $self->{'ID'}       = ($args{-ID} or 
			 $self->{ID} or 
			 "Unknown"); 
    $self->{'name'}     = ($args{-name} or 
			 $self->{name} or
			 "Unknown");
    $self->{'class'} = ($args{-class} or 
			 $self->{class} or
			 "Unknown");
    $self->{'strand'} = ($args{-strand} or 
			 $self->{strand} or
			 "+");
    $self->{'tags'} = $args{-tags} ? ((ref($args{-tags}) eq "HASH") ? $args{-tags} : {} ) :{};
    return $self;
}

=head2 ID

 Title   : ID
 Usage   : my $ID = $icm->ID()
           $pfm->ID('M00119');
 Function: Get/set on the ID of the pattern (unique in a DB or a set)
 Returns : pattern ID (a string)
 Args    : none for get, string for set

=cut

sub ID  {
    my ($self, $ID) = @_;
    $self->{'ID'} = $ID if $ID;
    return $self->{'ID'};
}


=head2 name

 Title   : name
 Usage   : my $name = $pwm->name()
           $pfm->name('PPARgamma');
 Function: Get/set on the name of the pattern
 Returns : pattern name (a string)
 Args    : none for get, string for set

=cut

sub name  {
    my ($self, $name) = @_;
    $self->{'name'} = $name if $name;
    return $self->{'name'};
}


=head2 class

 Title   : class
 Usage   : my $class = $pwm->class()
           $pfm->class('forkhead');
 Function: Get/set on the structural class of the pattern
 Returns : class name (a string)
 Args    : none for get, string for set

=cut


sub class  {
    my ($self, $class) = @_;
    $self->{'class'} = $class if $class;
    return $self->{'class'};
}

=head2 tag

 Title   : tag
 Usage   : my $acc = $pwm->tag('acc')
           $pfm->tag(source => "Gibbs");
 Function: Get/set on the structural class of the pattern
 Returns : tag value (a scalar/reference)
 Args    : tag name (string) for get,
	   tag name (string) and value (any scalar/reference) for set

=cut

sub tag  {
    my ($self, $tag, $value) = @_;
    return unless $tag;
    if (defined $value) {
	$self->{'tags'}->{$tag} = $value;
    }
    return $self->{'tags'}->{$tag};
}

=head2 all_tags

 Title   : all_tags
 Usage   : my %tag = $pfm->all_tags();
 Function: get a hash of all tags for a matrix
 Returns : a hash of all tag values keyed by tag name
 Args    : none

=cut

sub all_tags {
    return %{$_[0]->{'tags'}};
}


=head2 matrix

 Title   : matrix
 Usage   : my $matrix = $pwm->matrix();
	   $pwm->matrix( [ [12, 3, 0, 0, 4, 0],
			   [ 0, 0, 0,11, 7, 0],
			   [ 0, 9,12, 0, 0, 0],
			   [ 0, 0, 0, 1, 1,12]
			 ]);

 Function: get/set for the matrix data
 Returns : a reference to 2D array of integers(PFM) or floats (ICM, PWM)
 Args    : none for get;
	   a four line string, reference to 2D array, or a 2D piddle for set

=cut


sub matrix  {
    my ($self, $matrixdata) = @_;
    $self->set_matrix($matrixdata) if $matrixdata;
    my @list = $self->{'matrix'}->list;
    my @array;
    for (0..3)  {
	push @array, [splice(@list, 0, $self->length())];
    }
    return \@array;

}

=head2 pdl_matrix

 Title   : pdl_matrix
 Usage   : my $pdl = $pwm->pdl_matrix();
 Function: access the PDL matrix used to store the actual
	   matrix data directly
 Returns : a PDL object, aka a piddle
 Args    : none

=cut

sub pdl_matrix  {
    $_[0]->{'Matrix'};
}

sub set_matrix  {
    my ($self, $matrixdata) = @_;

    # The input matrix (specified as -array=> in the constructir call
    # can either be 
    #      * a 2D regular perl array with 4 rows,
    #      * a piddle (FIXME - check for 4 rows), or
    #      * a four-line string of numbers

    # print STDERR "MATRIX>>>".$matrixdata;
    if (ref($matrixdata) eq "ARRAY"         
	and ref($matrixdata->[0]) eq "ARRAY"
	and scalar(@{$matrixdata}) == 4)  
    {  
        # it is a perl array
	$self->{matrix} = pdl $matrixdata;
    }
    elsif (ref($matrixdata) eq "PDL")  
    {   
        # it's a piddle
	$self->{matrix} = $matrixdata;
    }
    elsif (!ref($matrixdata))
	   #and (scalar split "\n",$matrixdata) == 4)
    {
	# it's a string then
	$self->{matrix} = $self->_matrix_from_string($matrixdata);
    }
    else  {
	$self->throw("Wrong data type/format for -matrix.\n".
		      "Acceptable formats are Array of Arrays (4 rows),\n".
		      "PDL Array, (4 rows),\n".
		      "or plain string (4 lines).");
    }
    # $self->_set_min_max_score();
    return 1;

}

sub _matrix_from_string  {
    my ($self, $matrixstring) = @_;
    my @array = ();
    foreach ((split "\n", $matrixstring)[0..3])  {
	s/^\s+//;
	s/\s+$//;
	push @array, [split];
    }
    return pdl \@array;
}

sub _set_min_max_score  {
    my ($self) = @_;
    my $transpose = $self->{matrix}->xchg(0,1);
    $self->{min_score} = sum(minimum $transpose);
    $self->{max_score} = sum(maximum $transpose);
}

sub _load {
    my ($self, $field, $value) = @_;
    if (substr(ref($self->{db}),0,5) eq "DBI::")  {
	# database retrieval
    }
    elsif (-d $self->{dbh})  {
	# retrieval from .pwm files in a directory
	$self->_lookup_in_matrixlist($field, $value) 
	    or do {
		warn ("Matrix with $field=>$value not found.");
		return undef;
	    };
	my $ID = $self->{ID};
	my $DIR = $self->{dbh};
	$self->set_matrix(scalar `cat $DIR/$ID.pwm`); # FIXME - temporary
	
    }
    else  {
	$self->throw("-dbh is not a valid database handle or a directory."); 
    }
}


=head2 revcom

 Title   : revcom
 Usage   : my $revcom_pfm = $pfm->revcom();
 Function: create a matrix pattern object which is reverse complement
	    of the current one
 Returns : a TFBS::Matrix::* object of the same type as the one
	    the method acted upon
 Args    : none

=cut

sub revcom  {
    my ($self) = @_;
    my $revcom_matrix = 
	$self->new(-matrix => $self->{matrix}->slice('-1:0,-1:0'),
		   # the above line rotates the original matrix 180 deg,
		   -ID       => ($self->{ID} or ""),
		   -name     => ($self->{name} or ""),
		   -class => ($self->{class} or ""),
		   -strand   => ($self->{strand} and $self->{strand} eq "-") ? "+" : "-",
		   -tags     => ($self->{tags} or {})  );
    return $revcom_matrix;
}


=head2 rawprint

 Title   : rawprint
 Usage   : my $rawstring = $pfm->rawprint);
 Function: convert matrix data to a simple tab-separated format
 Returns : a four-line string of tab-separated integers or floats
 Args    : none

=cut


sub rawprint  {
    my $self = shift;
    my $pwmstring = sprintf ( $self->{matrix} );
    $pwmstring =~ s/\[|\]//g;                # lose []
    $pwmstring =~ s/\n /\n/g;                # lose leading spaces
    my @pwmlines = split("\n", $pwmstring); # f
    $pwmstring = join ("\n", @pwmlines[2..5])."\n";
    return $pwmstring;
}

=head2 prettyprint

 Title   : prettyprint
 Usage   : my $prettystring = $pfm->prettyprint();
 Function: convert matrix data to a human-readable string format
 Returns : a four-line string with nucleotides and aligned numbers
 Args    : none

=cut

sub prettyprint  {
    my $self = shift;
    my $pwmstring = sprintf ( $self->{matrix} );
    $pwmstring =~ s/\[|\]//g;                # lose []
    $pwmstring =~ s/\n /\n/g;                # lose leading spaces
    my @pwmlines = split("\n", $pwmstring); # 
    @pwmlines = ("A  [$pwmlines[2] ]",
		 "C  [$pwmlines[3] ]",
		 "G  [$pwmlines[4] ]",
		 "T  [$pwmlines[5] ]");
    $pwmstring = join ("\n", @pwmlines)."\n";
    return $pwmstring;
}

=head2 length

 Title   : length
 Usage   : my $pattern_length = $pfm->length;
 Function: gets the pattern length in nucleotides
	    (i.e. number of columns in the matrix)
 Returns : an integer
 Args    : none

=cut

sub length  {
    my $self = shift;
    return $self->{matrix}->getdim(0);
}


sub DESTROY  {
    # nothing
}



1;





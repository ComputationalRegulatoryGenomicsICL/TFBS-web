# TFBS module for TFBS::MatrixSet
#
# Copyright Boris Lenhard
#
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

TFBS::Matrix::Set - an agregate class representing a set of matrix patterns, containing methods for manipulating the set as a whole

=head1 SYNOPSIS

    # creation of a TFBS::MatrixSet object
    # let @list_of_matrix_objects be a list of TFBS::Matrix::* objects
 
    ###################################
    # Create a TFBS::MatrixSet object:
    
    my $matrixset = TFBS::MatrixSet->new(); # creates an empty set
    $matrixset->add_Matrix(@list_of_matrix_objects); #add matrix objects to set
    $matrixset->add_Matrix($matrixobj); # adds a single matrix object to set

    # or, same as above:
    
    my $matrixset = TFBS::MatrixSet->new(@list_of_matrix_objects, $matrixobj);

    ###################################
    #


=head1 DESCRIPTION

TFBS::MatrixSet is an aggregate class storing a set of TFBS::Matrix::* subclass objects, and providing methods form manipulating those sets as a whole. TFBS::MatrixSet objects are created <I>de novo<I> or returned by some database (TFBS::DB::*) retrieval methods.

=head1 FEEDBACK

Please send bug reports and other comments to the author.

=head1 AUTHOR - Boris Lenhard

Boris Lenhard E<lt>Boris.Lenhard@cgb.ki.seE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are preceded with an underscore.

=cut



# The code begins HERE:

package TFBS::MatrixSet;
use vars '@ISA';

use PDL;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;
use File::Temp qw/:POSIX/;

use TFBS::Matrix;
use TFBS::SiteSet;

use strict;

@ISA = qw(Bio::Root::RootI);


=head2 new

=cut


sub new  {
    my ($caller, @matrices) = shift;
    my $self = bless {}, ref($caller) || $caller;
    $self->add_matrix(@matrices) if @matrices;
    return $self;
}


sub add_matrix  {
    my ($self, @matrices) = @_;
    push @{$self->{matrix_list}}, @matrices;
}

sub add_matrix_set  {
    my ($self, @sets) = @_;
    foreach my $matrixset (@sets)  {
	push @{$self->{matrix_list}}, @{$matrixset->{matrix_list}};
    }
}

sub reset {
    my ($self) = @_;
    @{$self->{_iterator_list}} = @{$self->{matrix_list}};
}

sub sort_by_name  {
    my ($self) = @_;
    @{$self->{matrix_list}} = sort { uc($a->{name}) cmp uc ($b->{name}) }
                              @{$self->{matrix_list}};
    $self->reset();
}

sub next {
    my ($self) = @_;
    if (my $next_matrix = shift (@{$self->{_iterator_list}})) {
	return $next_matrix;
    }
    else  { 
	$self->reset; 
	return undef;
    }
}

sub search_seq  {
    my ($self, %args) = @_;
    $self->_search(%args);
}


sub _search  { 

    my ($self, %args) = @_;

    # DIRTY - stick tmp file name to seq object

    my $seqobj = $self->_to_seqobj(%args);
    ($seqobj->{_fastaFH}, $seqobj->{_fastafile}) = tmpnam(); 
    # we need $fastafile below

    my $outstream = Bio::SeqIO->new(-file=>">".$seqobj->{_fastafile}, -format=>"Fasta");
    $outstream->write_seq($seqobj);
    $outstream->close;
    
    # iterate through pwms
    my @PWMs;
    $self->reset();

    while (my $pwm = $self->next() ) {
	push @PWMs,$pwm;
    }
    
    # do the analysis
   
    my $hitlist = TFBS::SiteSet->new();
    
    foreach my $pwm (@PWMs)  {
	my $threshold = ($args{-threshold} or $pwm->{minscore});
	$hitlist->add_siteset($pwm->search_seq(-seqobj=>$seqobj, 
					    -threshold =>$threshold ));
    }
    delete $seqobj->{_fastaFH};
    unlink $seqobj->{_fastafile};
    delete $seqobj->{_fastafile};
    return $hitlist;
}

sub _csearch  { 

    my ($self, %args) = @_;
    my $PWM_SEARCH = '/home/httpd/cgi-bin/CONSITE/bin/pwm_searchPFF';

    # DIRTY - stick tmp file name to seq object

    my $seqobj = $self->_to_seqobj(%args);
    ($seqobj->{_fastaFH}, $seqobj->{_fastafile}) = tmpnam(); 
    # we need $fastafile below

    my $seqFH = Bio::SeqIO->newFh(-fh=>$seqobj->{_fastaFH}, -format=>"Fasta");
    print $seqFH $seqobj;
 
    
    # iterate through pwms
    my @PWMs;
    $self->reset();

    while (my $pwm = $self->next() ) {
	push @PWMs,$pwm;
    }
    
    # do the analysis
   
    my $hitlist = TFBS::SiteSet->new();
    
    foreach my $pwm (@PWMs)  {
	my $threshold = ($args{-threshold} or $pwm->{minscore});
	$hitlist->add_siteset($pwm->search_seq(-seqobj=>$seqobj, 
					    -threshold =>$threshold ));
    }
    delete $seqobj->{_fastaFH};
    delete $seqobj->{_fastafile};
    return $hitlist;

}    



sub _bsearch  {  
    my ($self,%args) = @_; #the rest of @_ goes to _to_seqob;
    my @PWMs;

    # prepare the sequence

    my $seqobj = $self->_to_seqobj(%args);
    $seqobj->{_pdl_matrix} = _seq_to_pdlmatrix($seqobj);    
    
    # prepare the PWMs

    $self->reset();

    while (my $pwm = $self->next() ) {
	push @PWMs,$pwm;
    }
    
    # do the analysis
   
    my $hitlist = TFBS::SiteSet->new();
    
    foreach my $pwm (@PWMs)  {
	my $threshold = ($args{-threshold} or $pwm->{minscore});
	$hitlist->add_siteset($pwm->bsearch(-seqobj=>$seqobj, 
					    -threshold =>$threshold ));
    }
    delete $seqobj->{_pdl_matrix};
    return $hitlist;
}

sub _seq_to_pdlmatrix  {

    # not OO - help function for search

    my $seqobj = shift;
    my $seqstring = uc($seqobj->seq());

    my @perlarray;
    foreach (qw(A C G T))  {
	my $seqtobits = $seqstring;
	eval "\$seqtobits =~ tr/$_/1/";  # curr. letter $_ to 1
	eval "\$seqtobits =~ tr/1/0/c";  # non-1s to 0
	push @perlarray, [split("", $seqtobits)];
    }
    return byte (\@perlarray);
}



sub _to_seqobj {
    my ($self, %args) = @_;

    my $seq;
    if ($args{-file})  {    # not a Bio::Seq
	return Bio::SeqIO->new(-file => $args{-file},
			     -format => 'fasta',
			     -moltype => 'dna')->next_seq();
    }
    elsif ($args{-seqstring} 
	   or $args{-seq}) 
    {   # I guess it's a string then
	return Bio::Seq->new(-seq  => ($args{-seqstring} or $args{-seq}),
			     -format => 'fasta',
			     -moltype => 'dna');
    }
    elsif ($args{'-seqobj'} and ref($args{'-seqobj'}) =~ /Bio\:\:Seq/) {
	# do nothing (maybe check later)
	return $args{'-seqobj'};
    }
    #elsif (ref($format) =~ /Bio\:\:Seq/ and !defined $seq)  {
	# if only one parameter passed and it's a Bio::Seq
	#return $format;
    #}
    else  {
	$self->throw ("Wrong parametes passed to search method: ".%args);
    }


    # CONTINUE HERE TOMORROW

}
	
sub add_Matrix {
    my ($self, @matrixlist)  = @_;
    foreach (@matrixlist)  {
	ref($_) =~ /TFBS::Matrix::/ 
	    or $self->throw("Attempted to add an element ".
			     "that is not a TFBS::Matrix object.");
	push @{$self->{matrix_list}},  $_;
	push @{$self->{_iterator_list}}, $_;

    }
    return 1;
}



1;


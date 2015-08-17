# TFBS module for TFBS::Matrix::PFM
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

               
TFBS::Matrix::PFM - class for raw position frequency matrix patterns


=head1 SYNOPSIS

=over 4

=item * creating a TFBS::Matrix::PFM object manually:


    my $matrixref = [ [ 12,  3,  0,  0,  4,  0 ],
		      [  0,  0,  0, 11,  7,  0 ],
		      [  0,  9, 12,  0,  0,  0 ],
		      [  0,  0,  0,  1,  1, 12 ]
		    ];	
    my $pfm = TFBS::Matrix::PFM->new(-matrix => $matrixref,
				     -name   => "MyProfile",
				     -ID     => "M0001"
				    );
    # or
 
    my $matrixstring =
        "12 3 0 0 4 0\n0 0 0 11 7 0\n0 9 12 0 0 0\n0 0 0 1 1 12";
 
    my $pfm = TFBS::Matrix::PFM->new(-matrixstring => $matrixstring,
				     -name   	   => "MyProfile",
				     -ID           => "M0001"
				    );
 
 
=item * retrieving a TFBS::Matix::PFM object from a database:

(See documentation of individual TFBS::DB::* modules to learn
how to connect to different types of pattern databases and 
retrieve TFBS::Matrix::* objects from them.)
    
    my $db_obj = TFBS::DB::JASPAR2->new
		    (-connect => ["dbi:mysql:JASPAR2:myhost",
				  "myusername", "mypassword"]);
    my $pfm = $db_obj->get_Matrix_by_ID("M0001", "PFM");
    # or
    my $pfm = $db_obj->get_Matrix_by_name("MyProfile", "PFM");


=item * retrieving list of individual TFBS::Matrix::PFM objects
from a TFBS::MatrixSet object

(See the L<TFBS::MatrixSet> to learn how to create 
objects for storage and manipulation of multiple matrices.)

    my @pfm_list = $matrixset->all_patterns(-sort_by=>"name");


=item * convert a raw frequency matrix to other matrix types:

    my $pwm = $pfm->to_PWM(); # convert to position weight matrix
    my $icm = $icm->to_ICM(); # convert to information con

=back        

=head1 DESCRIPTION

TFBS::Matrix::PFM is a class whose instances are objects representing
raw position frequency matrices (PFMs). A PFM is derived from N
nucleotide patterns of fixed size, e.g. the set of sequences

    AGGCCT
    AAGCCT
    AGGCAT
    AAGCCT
    AAGCCT
    AGGCAT
    AGGCCT
    AGGCAT
    AGGTTT
    AGGCAT
    AGGCCT
    AGGCCT


will give the matrix:

    A:[ 12  3  0  0  4  0 ]
    C:[  0  0  0 11  7  0 ]
    G:[  0  9 12  0  0  0 ]
    T:[  0  0  0  1  1 12 ]

which contains the count of each nucleotide at each position in the
sequence. (If you have a set of sequences as above and want to
create a TFBS::Matrix::PFM object out of them, have a look at
TFBS::PatternGen::SimplePFM module.)

PFMs are easily converted to other types of matrices, namely
information content matrices and position weight matrices. A
TFBS::Matrix::PFM object has the methods to_ICM and to_PWM which
do just that, returning a TFBS::Matrix::ICM and TFBS::Matrix::PWM
objects, respectively. 

=head1 FEEDBACK

Please send bug reports and other comments to the author.

=head1 AUTHOR - Boris Lenhard

Boris Lenhard E<lt>Boris.Lenhard@cgb.ki.seE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are preceded with an underscore.

=cut


# The code begins HERE:

package TFBS::Matrix::PFM;

use vars '@ISA';
use PDL;
use strict;
use Bio::Root::RootI;
use Bio::SeqIO;
use TFBS::Matrix;
use TFBS::Matrix::ICM;
use TFBS::Matrix::PWM;
use File::Temp qw/:POSIX/;
@ISA = qw(TFBS::Matrix Bio::Root::RootI);

#######################################################
# PUBLIC METHODS
#######################################################

=head2 new

 Title   : new
 Usage   : my $pfm = TFBS::Matrix::PFM->new(%args)
 Function: constructor for the TFBS::Matrix::PFM object
 Returns : a new TFBS::Matrix::PFM object
 Args    : # you must specify either one of the following three:
 
	   -matrix,      # reference to an array of arrays of integers
	      #or
	   -matrixstring,# a string containing four lines
	                 # of tab- or space-delimited integers
	      #or
	   -matrixfile,  # the name of a file containing four lines
	                 # of tab- or space-delimited integers
	   #######
 
           -name,        # string, OPTIONAL
           -ID,          # string, OPTIONAL
           -class,       # string, OPTIONAL
           -tags         # an array reference, OPTIONAL
Warnings  : Warns if the matrix provided has columns with different
            sums. Columns with different sums contradict the usual
	    origin of matrix data and, unless you are absolutely sure
	    that column sums _should_ be different, it would be wise to
	    check your matrices.

=cut

sub new  {
    my ($class, %args) = @_;
    my $matrix = TFBS::Matrix->new(%args, -matrixtype=>"PFM");
    my $self = bless $matrix, ref($class) || $class;
    $self->_check_column_sums();
    return $self;
}

=head2 column_sum

 Title   : column_sum
 Usage   : my $nr_sequences = $pfm->column_sum()
 Function: calculates the sum of elements of one column
	   (the first one by default) which normally equals the
           number of sequences used to derive the PFM. 
 Returns : the sum of elements of one column (an integer)
 Args    : columnn number (starting from 1), OPTIONAL - you DO NOT
           need to specify it unless you are dealing with a matrix

=cut

sub column_sum {
    my ($self, $column) = (@_,1);
    return $self->{matrix}->slice($column-1)->sum;
    
}

=head2 to_PWM

 Title   : to_PWM
 Usage   : my $pwm = $pfm->to_PWM()
 Function: converts a raw frequency matrix (a TFBS::Matrix::PFM object)
	   to position weight matrix. At present it assumes uniform
	   background distribution of nucleotide frequencies.
 Returns : a new TFBS::Matrix::PWM object
 Args    : none; in the future releases, it should be able to accept
	   a user defined background probability of the four
	   nucleotides

=cut

sub to_PWM  {
    my ($self) = @_;
    my $nseqs = $self->{'matrix'}->sum / $self->length;
    my $q_pdl = ($self->{'matrix'} +0.25*sqrt($nseqs))
		 / 
		($nseqs + sqrt($nseqs));
    my $pwm_pdl = log2(4*$q_pdl);

    my $PWM = TFBS::Matrix::PWM->new
	( (map {("-$_", $self->{$_}) } keys %$self),
	  -matrix    => $pwm_pdl
	);
    return $PWM;
    
}


=head2 to_ICM

 Title   : to_ICM
 Usage   : my $icm = $pfm->to_ICM()
 Function: converts a raw frequency matrix (a TFBS::Matrix::PFM object)
	   to information content matrix. At present it assumes uniform
	   background distribution of nucleotide frequencies.
 Returns : a new TFBS::Matrix::ICM object
 Args    : none; in the future releases, it should be able to accept
	   a user defined background probability of the four
	   nucleotides

=cut

sub to_ICM  {
    my ($self) = @_;
    my $p_pdl = $self->{'matrix'} / $self->{'matrix'}->xchg(0,1)->sumover;
    my $plog_pdl = $p_pdl*log2($p_pdl);
    $plog_pdl = $plog_pdl->badmask(0);
    my $D_pdl = 2 + $plog_pdl->xchg(0,1)->sumover;
    my $ic_pdl = $p_pdl * $D_pdl;

    my $ICM = TFBS::Matrix::ICM->new
	( (map {("-$_" => $self->{$_})} keys %$self),
	  -matrix    => $ic_pdl
	);
    return $ICM;

}


=head2 draw_logo

 Title   : draw_logo
 Usage   : my $gd_image = $pfm->draw_logo()
 Function: draws a sequence logo; this is a shorcut function
	   that internally converis $pfm to an information
	   content matrix (TFBS::Matrix::ICM) object, and calls its
	   draw_logo function; for details how to use this method,
	   see draw_logo entry in TFBS::Matrix::ICM documentation
 Returns : a GD image object (see documentation of GD module)
 Args    : many;
	   see draw_logo entry in TFBS::Matrix::ICM documentation

=cut

sub draw_logo {
    my ($self, %args) = @_;
    $self->to_ICM->draw_logo(%args);
}

=head2 name

=head2 ID

=head2 class

=head2 matrix

=head2 length

=head2 revcom

=head2 rawprint

=head2 prettyprint

The above methods are common to all matrix objects. Please consult
L<TFBS::Matrix> to find out how to use them.

=cut

###############################################
# PRIVATE METHODS
###############################################

sub _check_column_sums  {
    my ($self) = @_;
    my $pdl = $self->{matrix}->sever();
    my $rowsums = $pdl->xchg(0,1)->sumover();
    if ($rowsums->where($rowsums != $rowsums->slice(0))->getdim(0) > 0)  {
	$self->warn("PFM for ".$self->{ID}." has unequal column sums");
    }
}

sub DESTROY  {
    # does nothing
}

###############################################
# UTILITY FUNCTIONS
###############################################

sub log2 { log($_[0]) / log(2); }


1;

# TFBS module for TFBS::Matrix::ICM
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

TFBS::Matrix::ICM - class for information content matrices of nucleotide
patterns

=head1 SYNOPSIS

=over 4

=item * creating a TFBS::Matrix::ICM object manually:

    my $matrixref = [ [ 0.00, 0.30, 0.00, 0.00, 0.24, 0.00 ],
		      [ 0.00, 0.00, 0.00, 1.45, 0.42, 0.00 ],
		      [ 0.00, 0.89, 2.00, 0.00, 0.00, 0.00 ],
		      [ 0.00, 0.00, 0.00, 0.13, 0.06, 2.00 ]
		    ];	
    my $icm = TFBS::Matrix::ICM->new(-matrix => $matrixref,
				     -name   => "MyProfile",
				     -ID     => "M0001"
				    );
 
    # or
 
    my $matrixstring = <<ENDMATRIX
    2.00   0.30   0.00   0.00   0.24   0.00
    0.00   0.00   0.00   1.45   0.42   0.00
    0.00   0.89   2.00   0.00   0.00   0.00
    0.00   0.00   0.00   0.13   0.06   2.00
    ENDMATRIX
    ;
    my $icm = TFBS::Matrix::ICM->new(-matrixstring => $matrixstring,
				     -name   	   => "MyProfile",
				     -ID           => "M0001"
				    );


=item * retrieving a TFBS::Matix::ICM object from a database:

(See documentation of individual TFBS::DB::* modules to learn
how to connect to different types of pattern databases and retrieve
TFBS::Matrix::* objects from them.)
    
    my $db_obj = TFBS::DB::JASPAR2->new
		    (-connect => ["dbi:mysql:JASPAR2:myhost",
				  "myusername", "mypassword"]);
    my $pfm = $db_obj->get_Matrix_by_ID("M0001", "ICM");
    # or
    my $pfm = $db_obj->get_Matrix_by_name("MyProfile", "ICM");


=item * retrieving list of individual TFBS::Matrix::ICM objects
from a TFBS::MatrixSet object

(see decumentation of TFBS::MatrixSet to learn how to create 
objects for storage and manipulation of multiple matrices)

    my @icm_list = $matrixset->all_patterns(-sort_by=>"name");

* drawing a sequence logo
          
    $icm->draw_logo(-file=>"logo.png", 
		    -full_scale =>2.25,
		    -xsize=>500,
		    -ysize =>250, 
		    -graph_title=>"C/EBPalpha binding site logo", 
		    -x_title=>"position", 
		    -y_title=>"bits");


=head1 DESCRIPTION

TFBS::Matrix::ICM is a class whose instances are objects representing
position weight matrices (PFMs). An ICM is normally calculated from a
raw position frequency matrix (see L<TFBS::Matrix::PFM>
for the explanation of position frequency matrices). For example, given
the following position frequency matrix,

    A:[ 12     3     0     0     4     0  ]
    C:[  0     0     0    11     7     0  ]
    G:[  0     9    12     0     0     0  ]
    T:[  0     0     0     1     1    12  ]

the standard computational procedure is applied to convert it into the
following information content matrix:

    A:[2.00  0.30  0.00  0.00  0.24  0.00]
    C:[0.00  0.00  0.00  1.45  0.42  0.00]
    G:[0.00  0.89  2.00  0.00  0.00  0.00]
    T:[0.00  0.00  0.00  0.13  0.06  2.00]

which contains the "weights" associated with the occurence of each
nucleotide at the given position in a pattern.

A TFBS::Matrix::PWM object is equipped with methods to search nucleotide
sequences and pairwise alignments of nucleotide sequences with the
pattern they represent, and return a set of sites in nucleotide
sequence (a TFBS::SiteSet object for single sequence search, and a
TFBS::SitePairSet for the alignment search).

=head1 FEEDBACK

Please send bug reports and other comments to the author.

=head1 AUTHOR - Boris Lenhard

Boris Lenhard E<lt>Boris.Lenhard@cgb.ki.seE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are preceded with an underscore.

=cut

# The code starts HERE:

package TFBS::Matrix::ICM;

use vars '@ISA';
use PDL;
use strict;
use Bio::Root::RootI;
use Bio::SeqIO;
use TFBS::Matrix;
use GD;
use File::Temp qw/:POSIX/;
@ISA = qw(TFBS::Matrix Bio::Root::RootI);

#################################################################
# PUBLIC METHODS
#################################################################

=head2 new

 Title   : new
 Usage   : my $icm = TFBS::Matrix::ICM->new(%args)
 Function: constructor for the TFBS::Matrix::ICM object
 Returns : a new TFBS::Matrix::ICM object
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

=cut


sub new  {
    my ($class, %args) = @_;
    my $matrix = TFBS::Matrix->new(%args, -matrixtype=>"ICM");
    my $self = bless $matrix, ref($class) || $class;
    $self->_check_ic_validity();
    return $self;
}

=head2 to_PWM

 Title   : to_PWM
 Usage   : my $pwm = $icm->to_PWM()
 Function: converts an  information content matrix (a TFBS::Matrix::ICM object)
	   to position weight matrix. At present it assumes uniform
	   background distribution of nucleotide frequencies.
 Returns : a new TFBS::Matrix::PWM object
 Args    : none; in the future releases, it should be able to accept
	   a user defined background probability of the four
	   nucleotides

=cut


sub to_PWM  {
    my ($self) = @_;
    $self->throw ("Method to_PWM not yet implemented.");
}

=head2 draw_logo

 Title   : draw_logo
 Usage   : my $gdImageObj = $icm->draw_logo(%args)
 Function: Draws a "sequence logo", a graphical representation
	   of a possibly degenerate fixed-width nucleotide
	   sequence pattern, from the information content matrix
 Returns : a GD::Image object;
	   if you only need the image file you can ignore it
 Args    : -file,       # the name of the output PNG image file
		        # OPTIONAL: default none
	   -xsize       # width of the image in pixels
		        # OPTIONAL: default 600
	   -ysize       # height of the image in pixels
		        # OPTIONAL: default 5/8 of -x_size
	   -margin      # size of image margins in pixels
		        # OPTIONAL: default 15% of -y_size
	   -full_scale  # the maximum value on the y-axis, in bits
		        # OPTIONAL: default 2.25
	   -graph_title,# the graph title
			# OPTIONAL: default none
	   -x_title,    # x-axis title; OPTIONAL: default none
	   -y_title     # y-axis title; OPTIONAL: default none


=cut

sub draw_logo {
    my $self = shift;
    my %args = (-xsize      => 600,
		-full_scale => 2.25,
		-graph_title=> "",
		-x_title    => "",
		-y_title    => "",
		@_);
    # Other parameters that can be specified:
    #       -ysize -line_width -margin
    # do not have a fixed default value 
    #   - they are calculated from xsize if not specified

    my ($xsize,$FULL_SCALE, $x_title, $y_title)   
	= @args{qw(-xsize -full_scale -x_title y_title)} ;

    my $PER_PIXEL_LINE = 300;
    
    # calculate other parameters if not specified

    my $line_width = ($args{-line_width} or int ($xsize/$PER_PIXEL_LINE) or 1);
    my $ysize      = ($args{-ysize} or $xsize/1.6); 
    # remark (the line above): 1.6 is a standard screen x:y ratio
    my $margin     = ($args{-margin} or $ysize*0.15);

    my $image = GD::Image->new($xsize, $ysize);
    my $white = $image->colorAllocate(255,255,255);
    my $black = $image->colorAllocate(0,0,0);
    my $motif_size = $self->{matrix}->getdim(0);
    my $font = ((gdTinyFont, gdSmallFont, gdMediumBoldFont, 
		gdLargeFont, gdGiantFont)[int(($ysize-50)/100)]
	or gdGiantFont);
    my $title_font = ((gdSmallFont, gdMediumBoldFont, 
		gdLargeFont, gdGiantFont)[int(($ysize-50)/100)]
	or gdGiantFont);


    # WRITE LABELS AND TITLE

    # graph title
    $image->string($title_font,
		   $xsize/2-length($args{-graph_title})*$title_font->width()/2,
		   $margin/2 - $title_font->height()/2,
		   $args{-graph_title}, $black);
    
    # x_title
    $image->string($font,
		   $xsize/2-length($args{-x_title})*$font->width()/2,
		   $ysize-( $margin - $font->height()*0 - 5*$line_width)/2 - $font->height()/2*0,
		   $args{-x_title}, $black);
    # y_title
    $image->stringUp($font,
		   ($margin -$font->width()- 5*$line_width)/2 - $font->height()/2 ,
		   $ysize/2+length($args{'-y_title'})*$font->width()/2,
		   $args{'-y_title'}, $black);
    

    # DRAW AXES

    # vertical: (top left to bottom right)
    $image->filledRectangle($margin-$line_width, $margin-$line_width, 
			 $margin-1, $ysize-$margin+$line_width, 
			 $black);
    # horizontal: (ditto)
    $image->filledRectangle($margin-$line_width, $ysize-$margin+1, 
			 $xsize-$margin+$line_width,$ysize-$margin+$line_width,
			 $black);

    # DRAW VERTICAL TICKS AND LABELS

    # vertical axis (IC 1 and 2) 
    my $ic_1 = ($ysize - 2* $margin) / $FULL_SCALE;
    foreach my $i (1..$FULL_SCALE)  {
	$image->filledRectangle($margin-3*$line_width, 
			     $ysize-$margin - $i*$ic_1, 
			     $margin-1, 
			     $ysize-$margin+$line_width - $i*$ic_1, 
			     $black);
	$image->string($font, 
		       $margin-5*$line_width - $font->width,
		       $ysize - $margin - $i*$ic_1 - $font->height()/2,
		       $i,
		       $black);
    }
    
    # DRAW HORIZONTAL TICKS AND LABELS, AND THE LOGO ITSELF 

    # define function refs as hash elements
    my %draw_letter = ( A => \&draw_A,
			C => \&draw_C,
			G => \&draw_G,
			T => \&draw_T );

    my $horiz_step = ($xsize -2*$margin) / $motif_size;
    foreach my $i (0..$motif_size)  {
	
	$image->filledRectangle($margin + $i*$horiz_step, 
			     $ysize-$margin+1, 
			     $margin + $i*$horiz_step+ $line_width, 
			     $ysize-$margin+3*$line_width, 
			     $black);
	last if $i==$motif_size;

	# get the $i-th column of matrix
	my %ic; 
	($ic{A}, $ic{C}, $ic{G}, $ic{T}) = list $self->{matrix}->slice($i);

	# sort nucleotides by increasing information content
	my @draw_order = sort {$ic{$a}<=>$ic{$b}} qw(A C G T);

	# draw logo column
	my $xlettersize = $horiz_step /1.1;
	my $ybottom = $ysize - $margin;
	foreach my $base (@draw_order)  {
	    my $ylettersize = int($ic{$base}*$ic_1 +0.5);
	    next if $ylettersize ==0;

	    # draw letter
	    $draw_letter{$base}->($image,
				  $margin + $i*$horiz_step,
				  $ybottom - $ylettersize,
			  $xlettersize, $ylettersize, $white);
	    $ybottom = $ybottom - $ylettersize-1;
	}	    
	    	
	$image->string($font,
		       $margin + ($i+0.5)*$horiz_step - $font->width()/2,
		       $ysize - $margin +5*$line_width,
		       $i+1,
		       $black);
    }

    # print $args{-file};
    if  ($args{-file}) {  
	open (PNGFILE, ">".$args{-file})
	    or $self->throw("Could not write to ".$args{-file});
        print PNGFILE $image->png;
	close PNGFILE;
    }
    return $image;
}

sub total_ic  {
    return $_[0]->{matrix}->sum();
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


#################################################################
# INTERNAL METHODS
#################################################################


sub _check_ic_validity  {
    my ($self) = @_;
    # to do
}

sub DESTROY  {
    # nothing
}


#################################################################
# UTILITY FUNCTIONS
#################################################################


# letter drawing routines

sub draw_A {
    
    my ($im, $x, $y, $xsize, $ysize, $white) = @_;
    my $green = $im->colorAllocate(0,255,0);
    my $outPoly = GD::Polygon->new();
    $outPoly->addPt($x, $y+$ysize);
    $outPoly->addPt($x+$xsize*.42, $y);
    $outPoly->addPt($x+$xsize*.58, $y);
    $outPoly->addPt($x+$xsize, $y+$ysize);
    $outPoly->addPt($x+0.85*$xsize, $y+$ysize);
    $outPoly->addPt($x+0.725*$xsize, $y+0.75*$ysize);
    $outPoly->addPt($x+0.275*$xsize, $y+0.75*$ysize);
    $outPoly->addPt($x+0.15*$xsize, $y+$ysize);
    $im->filledPolygon($outPoly, $green);
    if ($ysize>8)  {
	my $inPoly = GD::Polygon->new();
	$inPoly->addPt($x+$xsize*.5, $y+0.2*$ysize);
	$inPoly->addPt($x+$xsize*.34, $y+0.6*$ysize-1);
	$inPoly->addPt($x+$xsize*.64, $y+0.6*$ysize-1);
	$im->filledPolygon($inPoly, $white);
    }
    return 1;
}
    
sub draw_C  {
    my ($im, $x, $y, $xsize, $ysize, $white) = @_;
    my $blue = $im->colorAllocate(0,0,255);
    $im->arc($x+$xsize*0.54, $y+$ysize/2,1.08*$xsize,$ysize,0,360,$blue);
    $im->fill($x+$xsize/2, $y+$ysize/2, $blue);
    if ($ysize>12) {
	$im->arc($x+$xsize*0.53, $y+$ysize/2, 
		 0.75*$xsize, (0.725-0.725/$ysize)*$ysize,
		 0,360,$white);
	$im->fill($x+$xsize/2, $y+$ysize/2, $white);
	$im->filledRectangle($x+$xsize/2, $y+$ysize/4+1, 
			     $x+$xsize*1.1, $y+(3*$ysize/4)-1,
			     $white);
    }
    elsif ($ysize>3)  {
	$im->arc($x+$xsize*0.53, $y+$ysize/2, 
		 (0.75-0.75/$ysize)*$xsize, (0.725-0.725/$ysize)*$ysize,
		 0,360,$white);
	$im->fill($x+$xsize/2, $y+$ysize/2, $white);
	$im->filledRectangle($x+$xsize*0.25, $y+$ysize/2, 
			     $x+$xsize*1.1, $y+$ysize/2,
			     $white);

    }
   return 1;
}

sub draw_G  {
    my ($im, $x, $y, $xsize, $ysize, $white) = @_;
    my $yellow = $im->colorAllocate(200,200,0);
    $im->arc($x+$xsize*0.54, $y+$ysize/2,1.08*$xsize,$ysize,0,360,$yellow);
    $im->fill($x+$xsize/2, $y+$ysize/2, $yellow);
    if ($ysize>20) {
	$im->arc($x+$xsize*0.53, $y+$ysize/2, 
		 0.75*$xsize, (0.725-0.725/$ysize)*$ysize,
		 0,360,$white);
	$im->fill($x+$xsize/2, $y+$ysize/2, $white);
	$im->filledRectangle($x+$xsize/2, $y+$ysize/4+1, 
			     $x+$xsize*1.1, $y+$ysize/2-1,
			     $white);
    }
    elsif($ysize>3)  {
	$im->arc($x+$xsize*0.53, $y+$ysize/2, 
		 (0.75-0.75/$ysize)*$xsize, (0.725-0.725/$ysize)*$ysize,
		 0,360,$white);
	$im->fill($x+$xsize/2, $y+$ysize/2, $white);
	$im->filledRectangle($x+$xsize*0.25, $y+$ysize/2, 
			     $x+$xsize*1.1, $y+$ysize/2,
			     $white);

    }
    $im->filledRectangle($x+0.85*$xsize, $y+$ysize/2,
			 $x+$xsize,$y+(3*$ysize/4)-1,
			  $yellow);
    $im->filledRectangle($x+0.6*$xsize, $y+$ysize/2,
			 $x+$xsize,$y+(5*$ysize/8)-1,
			  $yellow);
   return 1;
}
    
sub draw_T {
    
    my ($im, $x, $y, $xsize, $ysize, $white) = @_;
    my $red = $im->colorAllocate(255,0,0);
    $im->filledRectangle($x, $y, $x+$xsize, $y+0.16*$ysize, $red);
    $im->filledRectangle($x+0.42*$xsize, $y, $x+0.58*$xsize, $y+$ysize, $red);
    return 1;
}
  

1;










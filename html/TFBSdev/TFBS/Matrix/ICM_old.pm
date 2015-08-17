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

sub new  {
    my ($class, %args) = @_;
    my $matrix = TFBS::Matrix->new(%args, -matrixtype=>"ICM");
    my $self = bless $matrix, ref($class) || $class;
    $self->_check_ic_validity();
    return $self;
}

sub _check_ic_validity  {
    my ($self) = @_;
    # to do
}

sub to_PWM  {
    my ($self) = @_;
    $self->throw ("Method to_PWM not yet implemented.");
}

sub draw_logo {
    my $self = shift;
    my %args = (-xsize      => 600,
		-full_scale => 2.25,
		-graph_title=> "",
		-x_title    => "",
		-y_title    => "",
		-letter_images => {A=>'/home/httpd/cgi-bin/CONSITE/A_univers.png',
				   C=>'/home/httpd/cgi-bin/CONSITE/C_univers.png',
				   G=>'/home/httpd/cgi-bin/CONSITE/G_univers.png',
				   T=>'/home/httpd/cgi-bin/CONSITE/T_univers.png'}, 
 
		@_);
    # Other parameters that can be specified:
    #       -ysize -line_width -margin
    # do not have a fixed default value 
    #   - they are calculated from xsize if not specified

    unless ($args{-file}) { $self->throw("No filename provided.") ;}

    my ($xsize,$FULL_SCALE, $x_title, $y_title)   
	= @args{qw(-xsize -full_scale -x_title y_title)} ;
    my %letter_file = %{$args{-letter_images}};
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
		   $ysize/2+length($args{-y_title})*$font->width()/2,
		   $args{-y_title}, $black);
    

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

    # load font images
    my %letter_image;
    foreach my $base(keys %letter_file) {
	$letter_image{$base} = GD::Image->newFromPng($letter_file{$base});
	$letter_image{$base}->transparent($white);
    }


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
	# print "\nPosition $i\n";
	# draw logo column
	my $xlettersize = $horiz_step;
	my $ybottom = $ysize - $margin +1;
	foreach my $base (@draw_order)  {
	    my $ylettersize = int($ic{$base}*$ic_1 +0.5);
	    next if $ylettersize ==0;

	    # draw letter

	    my ($LETTER_XSIZE, $LETTER_YSIZE) = 
		$letter_image{$base}->getBounds();
  	    $image->copyResized($letter_image{$base}, 
  				## destination x, y :
  				$margin + $i*$horiz_step,
  				$ybottom - $ylettersize,
  				## source x, y :
  				0, 0,
  				## destination width, height:
  				$xlettersize, $ylettersize,
  				## source width, height :
  				$LETTER_XSIZE, $LETTER_YSIZE);

	    # fix for squashed letters:

	    # print "ylettersize $ylettersize\n";
	    if ($ylettersize < 10)  {
		# middle fix
  		$image->copyResized($letter_image{$base}, 
  				    ## destination x, y :
  				    $margin + $i*$horiz_step,
  				    $ybottom - $ylettersize +1,
  				    ## source x, y :
  				    0, $LETTER_YSIZE*1/4,
  				    ## destination width, height:
  				    $xlettersize, $ylettersize-2,
  				    ## source width, height :
  				    $LETTER_XSIZE, $LETTER_YSIZE*3/4);
		# bottom fix - bottom first because of T
  		$image->copyResized($letter_image{$base}, 
  				    ## destination x, y :
  				    $margin + $i*$horiz_step,
  				    $ybottom - 1,
  				    ## source x, y :
  				    0, $LETTER_YSIZE*13/15,
  				    ## destination width, height:
  				    $xlettersize,  1,
  				    ## source width, height :
  				    $LETTER_XSIZE, $LETTER_YSIZE/15);
		# top fix
  		$image->copyResized($letter_image{$base}, 
  				    ## destination x, y :
  				    $margin + $i*$horiz_step,
  				    $ybottom - $ylettersize,
  				    ## source x, y :
  				    0, $LETTER_YSIZE*0.4/10,
  				    ## destination width, height:
  				    $xlettersize, 1,
  				    ## source width, height :
  				    $LETTER_XSIZE, $LETTER_YSIZE/15);
	    }
	    $ybottom = $ybottom - $ylettersize;
	}	    
	    	
	$image->string($font,
		       $margin + ($i+0.5)*$horiz_step - $font->width()/2,
		       $ysize - $margin +5*$line_width,
		       $i+1,
		       $black);
    }

    # print $args{-file};
    open (PNGFILE, ">".$args{-file}) or die;
    print PNGFILE $image->png;
    close PNGFILE;
    
}

sub total_ic  {
    return $_[0]->{matrix}->sum();
}

sub DESTROY  {
    # nothing
}


1;










package TFBS::Ext::pwmsearch;

require 5.005_62;
use strict;
use warnings;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Bio::SeqIO;
use File::Temp qw (:POSIX);

require Exporter;
require DynaLoader;

our @ISA = qw(Exporter DynaLoader);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use TFBS::Ext::pwmsearch ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
%EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

@EXPORT = qw(
	
);
$VERSION = '0.2';

bootstrap TFBS::Ext::pwmsearch $VERSION;

# Preloaded methods go here.

sub pwmsearch {
    my ($matrixobj, $seqobj, $threshold) = @_;

    my $matrixfile = tmpnam();
    open (MATRIX, ">$matrixfile") or die ("Error opening temporary file.");
    print MATRIX $matrixobj->rawprint();
    close MATRIX;

    my $outfile = tmpnam();

    my $seqfile;
    if ($seqobj->{_fastafile})  {
	$seqfile = $seqobj->{_fastafile};
    }
    else  {
	$seqfile = tmpnam();
	my $outstream = Bio::SeqIO->new(-file=>">$seqfile", -format=>"fasta");
	$outstream->write_seq($seqobj);
	$outstream->close();
    }
    # calculate threshold

    if ($threshold)  {
	if ($threshold =~ /(.+)%/)  { 
	    # percentage
	    $threshold = $matrixobj->{min_score} +
		($matrixobj->{max_score} - $matrixobj->{min_score})* $1/100;
	}
	else  {
	    # absolute value
	    # $threshold = $args{-threshold};
	}
    }
    else {
	# no threshold given
	$threshold = $matrixobj->{min_score} -1;
    }
    
    search_xs($matrixfile, $seqfile, 
	    $threshold, $matrixobj->name()."", 
	    $matrixobj->{'class'}."", $outfile);
    
    unlink $seqfile unless $seqobj->{'_fastafile'}; 
    unlink $matrixfile;  

    my $hitlist = TFBS::SiteSet->new();
    my ($TFname, $TFclass) = ($matrixobj->{name}, $matrixobj->{class});


    open (OUTFILE, $outfile) 
	or die("Could not read temporary outfile"); 
    while (my $line = <OUTFILE>)  {
	chomp $line;
	$line =~ s/^\s+//;
	my ($seq_id, $factor, $class, $strand, $score, $pos, $siteseq) =
	    (split /\s+/, $line)[0, 2, 3, 4, 5, 7, 9];
	my $num_strand = ($strand eq "-")? "-1" : "1";
	my $site = TFBS::Site->new ( -seqname => $seqobj->display_id()."",
				     -seqobj  => $seqobj,
				     -strand  => $num_strand."",
				     -pattern => $matrixobj,
				     -siteseq => $siteseq."",
				     -score   => $score."",
				     -start   => $pos,
				     -end     => $pos + length($siteseq) -1
				     );
	$hitlist->add_site($site);
    }
    close OUTFILE;
    unlink $outfile;
    return $hitlist;
}    


1;
__END__

=head1 NAME

TFBS::Ext::pwmsearch - Perl extension for scanning a DNA sequence object with a position weight matrix

=head1 SYNOPSIS

  use TFBS::Ext::pwmsearch;
  pwmsearch

=head1 DESCRIPTION

Stub documentation for TFBS::Ext::pwmsearch, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.


=head1 AUTHOR

A. U. Thor, a.u.thor@a.galaxy.far.far.away

=head1 SEE ALSO

perl(1).

=cut

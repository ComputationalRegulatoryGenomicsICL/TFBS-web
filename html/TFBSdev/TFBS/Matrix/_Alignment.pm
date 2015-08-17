package TFBS::Matrix::_Alignment;

use vars qw(@ISA $AUTOLOAD);

use TFBS::SitePair;
use TFBS::SitePairSet;
use Bio::Root::RootI;
use Bio::Seq;
use Bio::SimpleAlign;
use IO::String;
use PDL;

use strict;

@ISA =('Bio::Root::RootI');


# CONSTANTS

use constant DEFAULT_WINDOW => 50;
use constant DEFAULT_CUTOFF => 70;
use constant DEFAULT_THRESHOLD => "80%";


sub new  {
    
    # this is ugly; OK, OK, I'll rewrite it as soon as I can  
    my ($caller, %args) = @_;
    my $self = bless {}, ref $caller || $caller;
    $self->window($args{-window} or DEFAULT_WINDOW);
    $self->_parse_alignment(%args);
    $self->seq1length(length(_strip_gaps($self->alignseq1())));
    $self->seq2length(length(_strip_gaps($self->alignseq2())));
    $self->conservation1($self->_calculate_conservation($self->window(),1));
    $self->conservation2($self->_calculate_conservation($self->window(),2));
    $self->_do_sitesearch
	(($args{-pattern_set} or $self->throw("No -matrixset parameter")), 
	 ($args{-threshold} or DEFAULT_THRESHOLD), 
	 ($self->cutoff($args{-cutoff} or DEFAULT_CUTOFF)));

    # $self->_set_start_end(%args);  # Maybe later... 

    return $self;
    
}



sub DESTROY {
    # empty
}

sub _parse_alignment {
    my ($self, %args) = @_;
    my ($seq1, $seq2, $start);
    my $alignstring;

    if (defined $args{'-alignstring'})  {
	$alignstring = $args{'-alignstring'};
    }
    elsif (defined $args{'-file'})  {
	$alignstring = _alignfile_to_string($args{'-file'});
    }
    elsif  (defined $args{-alignobj})  {
	$alignstring = _alignobj_to_string($args{'alignobj'});
    }
    else  {
	$self->throw("No -alignstring, -file or -alignobj passed.");
    }
	
    
    my @match;
    my @alnlines = split("\n", $alignstring);
    shift @alnlines;shift @alnlines;shift @alnlines; # drop header
    
    while ($_=shift @alnlines) {
	$start=1;
	# print $_; 
	my ($label1, $string1) = split; # (/\s+/($_, 21,60);
	$self->seq1name($label1) unless $self->seq1name();
	my ($label2, $string2) = split /\s+/, shift(@alnlines); 
	$self->seq2name($label2) unless $self->seq2name();
                     #substr($string2,21,60);
	shift @alnlines; shift @alnlines; # skip asterisk line and blank line
	$seq1 .= $string1;
	$seq2 .= $string2;
    }

    $self->alignseq1($seq1);
    $self->alignseq2($seq2);
    my @seq1 = ("-", split('', $seq1) );
    my @seq2 = ("-", split('', $seq2) );
    $self->{alignseq1array} = [@seq1];
    $self->{alignseq2array} = [@seq2];
    
    my (@seq1index, @seq2index);
    my ($i1, $i2) = (0, 0);
    for my $pos (0..$#seq1) {
	my ($s1, $s2) = (0, 0);
	$seq1[$pos] ne "-" and  $s1 = ++$i1;
	$seq2[$pos] ne "-" and  $s2 = ++$i2;
	push @seq1index, $s1;
	push @seq2index, $s2;
    }

    $self->pdlindex( pdl [ [list sequence($#seq1+1)], 
			   [@seq1index], 
			   [@seq2index], 
			   [list zeroes ($#seq1+1)] ]) ;

#      # mark exons in row 3: only if sequence 1 contains exon objects
#      my $seqobj =$self->{job}->seq1();
#      my @features = $seqobj->top_SeqFeatures;
#      for my $exon (@features)  {
#  	next unless ref($exon) =~ /::Exon/;
#  	my $exonstart = $self->pdlindex($exon->start(), 1=>0);
#  	my $exonend   = $self->pdlindex($exon->end(),   1=>0);

#  	my $sl = $self->pdlindex->slice("$exonstart:$exonend,3");
#  	$sl++;
#      }

#      # mark ORF: only if cDNA contains generic feature with an 'orf' tag
    
#      my @cdnafeatures = ();
#      if ($self->{job}->seq3())  {
#  	@cdnafeatures = $self->{job}->seq3()->top_SeqFeatures;
#      }
#      my ($cdnaexonstart, $cdnaexonend) = (0,0);
#      my $orf;
#      for my $feat (@cdnafeatures)  {
#  	$orf = $feat 
#  	    if (ref($feat) eq "Bio::SeqFeature::Generic"
#  		and $feat->primary_tag eq 'orf');
#  	next unless ref($feat) =~ /::Exon/;
#  	$cdnaexonstart = $feat->start() if $cdnaexonstart > $feat->start();
#  	$cdnaexonend   = $feat->end()   if $cdnaexonend   < $feat->end();
#      }
#      if ($orf)  {
#  	my ($orf_start, $orf_end, $orf_strand) = 
#  	    ($orf->start(), $orf->end(), $orf->strand());
#  	my $cdnaPdl = 
#  	    $self->pdlindex->slice(':,3')->where
#  		($self->pdlindex->slice(':,3')==1);
#  	my $orfPdl;

#  	$orf_end = $cdnaexonend if $orf_end>$cdnaexonend;
#  	$orf_start = $cdnaexonstart if $orf_start < $cdnaexonstart;

#  	if ($orf_strand eq "1") {
#  	    $orf_end = $cdnaexonend if $orf_end>$cdnaexonend;
#  	    $orfPdl = $cdnaPdl->slice
#  		(($orf_start-$cdnaexonstart).":".($orf_end-$cdnaexonstart-1));
#  	}
#  	else  {
#  	    $orfPdl = $cdnaPdl->slice
#  		(($self->{job}->seq3->length() - $orf_end + $cdnaexonstart).
#  		 ":".
#  		 ($self->{job}->seq3->length()-$orf_start+$cdnaexonstart-1));
#  	}
#  	$orfPdl++; # set it to 2
#      }
	
    return 1;

}

sub pdlindex {
    my ($self, $input, $p1, $p2) = @_ ;
    # print ("PARAMS ", join(":", @_), "\n");
    if (ref($input) eq "PDL")  {
	$self->{pdlindex} = $input;
    }

    unless (defined $p2)  {
	return $self->{pdlindex};
    }
    else {
	my @results = list 
	    $self->{pdlindex}->xchg(0,1)->slice($p2)->where
		($self->{pdlindex}->xchg(0,1)->slice($p1)==$input);
	wantarray ? return @results : return $results[0];
    }
}

sub lower_pdlindex {
    my ($self, $input, $p1, $p2) = @_;
    unless (defined $p2)  {
	$self->throw("Wrong number of parameters passed to lower_pdlindex");
    }
    my $result;
    my $i = $input;

    until ($result = $self->pdlindex($i, $p1 => $p2))  { 
	$i--;

	last if $i==0;
    }
    return $result or 1;
}

sub higher_pdlindex {
    my ($self, $input, $p1, $p2) = @_;
    unless (defined $p2)  {
	$self->throw("Wrong number of parameters passed to lower_pdlindex");
    }
    my $result;
    my $i = $input;
    until ($result = $self->pdlindex($i, $p1 => $p2))  { 
	$i++; 
	last unless ($self->pdlindex($i, $p1=>0) > 0);
    }
    return $result;
}

			   
sub _calculate_conservation  {
    my ($self, $WINDOW, $which) = @_;
    my (@seq1, @seq2);
    if ($which==2)  {
	@seq1 = @{$self->{alignseq2array}};
	@seq2 = @{$self->{alignseq1array}};
    }
    else  {
	@seq1 = @{$self->{alignseq1array}};
	@seq2 = @{$self->{alignseq2array}};
	$which=1;
    }
    
    my @CONSERVATION;
    my @match;

    while ($seq1[0] eq "-")  {
	shift @seq1;
	shift @seq2;
    }

    for my $i (0..$#seq1) {
  	push (@match,( uc($seq1[$i]) eq uc($seq2[$i]) ? 1:0)) 
  	    unless ($seq1[$i] eq "-" or $seq1[$i] eq ".");
    }
    my @graph=($match[0]);
    for my $i (1..($#match+$WINDOW/2))  {
  	$graph[$i] = ($graph[$i-1] or 0) 
  	           + ($i>$#match ? 0: $match[$i]) 
  		   - ($i<$WINDOW ? 0: $match[$i-$WINDOW]);
    }

    # at this point, the graph values are shifted $WINDOW/2 to the right
    # i.e. the score at a certain position is the score of the window 
    # UPSTREAM of it: To fix it, we shoud discard the first $WINDOW/2 scores:
    #$self->conservation1 ([]);
    foreach my $pos (@graph[int($WINDOW/2)..$#graph])  {
	push @CONSERVATION, 100*$pos/$WINDOW;
    }

    return [@CONSERVATION];
    
}


sub _strip_gaps {
    # a utility function
    my $seq = shift;
    $seq =~ s/\-|\.//g;
    return $seq;
}


sub _do_sitesearch  {
    my ($self, $MATRIXSET, $THRESHOLD, $CUTOFF) = @_;
    
    my $seqobj1 = Bio::Seq->new(-seq=>_strip_gaps($self->alignseq1()),
				-id => "Seq1");
    my $siteset1 = 
	$MATRIXSET->search_seq(-seqobj => $seqobj1,
			    -threshold => $THRESHOLD);
    my $siteset1_itr = $siteset1->Iterator(-sort_by => "start");

    my $seqobj2 = Bio::Seq->new(-seq=>_strip_gaps($self->alignseq2()),
				-id => "Seq2");
    my $siteset2 =
	$MATRIXSET->search_seq(-seqobj => $seqobj2,
			    -threshold => $THRESHOLD);
    my $siteset2_itr = $siteset2->Iterator(-sort_by => "start");
    
    my $site1 = $siteset1_itr->next();
    my $site2 = $siteset2_itr->next();

    $self->site_pair_set(TFBS::SitePairSet->new()); 

    while (defined $site1 and defined $site2) {
	my $pos1_in_aln = $self->pdlindex($site1->start(), 1=>0);
	my $pos2_in_aln = $self->pdlindex($site2->start(), 2=>0);
	my $cmp = (($pos1_in_aln <=> $pos2_in_aln) 
		   or 
		   ($site1->pattern->name() cmp $site2->pattern->name()));

	if ($cmp==0) { ### match
	    if (# threshold test:
		$self->conservation1->[$site1->start()]
		>=
		$self->cutoff()
		#and 
		# exclude ORF test
		#(!($self->{job}->exclude_orf() 
		#   and $self->pdlindex($site1->start(), 1=>3) == 2
		#   )
		#)
		)
	    {
		#$self->fsiteset1->add_site($site1);
		#$self->fsiteset2->add_site($site2);
		my $site_pair = TFBS::SitePair->new($site1, $site2);
		$self->site_pair_set->add_site_pair($site_pair);
	    }
	    $site1 = $siteset1_itr->next();
	    $site2 = $siteset2_itr->next();
	}
	elsif ($cmp<0)  { ### $siteset1 is behind
	    $site1 = $siteset1_itr->next();
	}
	elsif ($cmp>0)  { ### $siteset2 is behind
	    $site2 = $siteset2_itr->next();
	}	    
    }    
}



sub _calculate_cutoff  {
    my ($self) = @_;
    my $ile = 0.9;
    my @conservation_array = sort {$a <=> $b} @{$self->conservation1()};

    my $perc_90 = $conservation_array[int($ile*scalar(@conservation_array))];
    return $perc_90;
}


sub _alignfile_to_string  {
    # a utility function
    my $alignfile = shift;
    if ($alignfile =~ /\.msf$/i) {
	my $alignobj = Bio::SimpleAlign->new();
	$alignobj->read_MSF($alignfile);
        return _alignobj_to_string($alignobj);
    }
    else  { #assumed clustalw - no AlignIO import yet
	local $/ = undef;
	open FILE, $alignfile
	    or die("Could not read alignfile $alignfile, stopped");
	my $alignstring = <FILE>;
	return $alignstring;
    }	
 }

sub _alignobj_to_string  {
    # a utility function
    my $alignobj = shift;
    my $alignstring;
    my $io = IO::String->new($alignstring);
    $alignobj->write_clustalw($io);
    return $alignstring;
}

# uglier than AUTOLOAD, but faster - a quick fix to get rid of Class::MethodMaker


sub cutoff          { $_[0]->{'cutoff'}        = $_[1] if exists $_[1]; $_[0]->{'cutoff'};       }
sub window          { $_[0]->{'window '}       = $_[1] if exists $_[1]; $_[0]->{'window '};      }
sub alignseq1       { $_[0]->{'alignseq1'}     = $_[1] if exists $_[1]; $_[0]->{'alignseq1'};    }
sub alignseq2       { $_[0]->{'alignseq2'}     = $_[1] if exists $_[1]; $_[0]->{'alignseq2'};    }
sub site_pair_set   { $_[0]->{'site_pair_set'} = $_[1] if exists $_[1]; $_[0]->{'site_pair_set'};}
sub seq1name        { $_[0]->{'seq1name'}      = $_[1] if exists $_[1]; $_[0]->{'seq1name'};     }
sub seq2name        { $_[0]->{'seq2name'}      = $_[1] if exists $_[1]; $_[0]->{'seq2name'};     }
sub seq1length      { $_[0]->{'seq1length'}    = $_[1] if exists $_[1]; $_[0]->{'seq1length'};   }
sub seq2length      { $_[0]->{'seq2length'}    = $_[1] if exists $_[1]; $_[0]->{'seq2length'};   }
sub conservation1   { $_[0]->{'conservation1'} = $_[1] if exists $_[1]; $_[0]->{'conservation1'};}
sub conservation2   { $_[0]->{'conservation2'} = $_[1] if exists $_[1]; $_[0]->{'conservation2'};}
sub exclude_orf     { $_[0]->{'exclude_orf'}   = $_[1] if exists $_[1]; $_[0]->{'exclude_orf'};  }
sub start_at        { $_[0]->{'start_at'}      = $_[1] if exists $_[1]; $_[0]->{'start_at'};     }
sub end_at          { $_[0]->{'end_at'}        = $_[1] if exists $_[1]; $_[0]->{'end_at'};     }
  



1;



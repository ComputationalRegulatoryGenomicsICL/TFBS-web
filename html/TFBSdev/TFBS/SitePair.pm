package TFBS::SitePair;

use vars qw(@ISA);
use strict;

use Bio::SeqFeature::FeaturePair;
@ISA = qw(Bio::SeqFeature::FeaturePair);

# 'new' is inherited 

sub pattern  {
    $_[0]->feature1->pattern();
}

sub site1  {
    $_[0]->feature1();
}

sub site2  {
    $_[0]->feature2();
}

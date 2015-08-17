package TFBS::SiteSet;

use vars '@ISA';
use TFBS::Site;
use strict;
@ISA = qw(Bio::Root::RootI);

sub new  {
    my ($class, @data) = @_;
    my $self = bless {}, ref($class) || $class;
    $self->{site_list} = $self->{_iterator_list} = [];
    @data = @{$class->{site_list}} if !@data && ref($class);
    $self->add_site(@data);
    return $self;
}

sub add_site {
    my ($self, @PWMlist)  = @_;
    foreach (@PWMlist)  {
	ref($_) =~ /TFBS::Site/ 
	    or $self->throw("Attempted to add an element ".
			     "that is not a TFBS::Site object.");
	push @{$self->{site_list}},  $_;
	push @{$self->{_iterator_list}}, $_;

    }
    return 1;
}

sub add_siteset {
    my ($self, @sitesets) = @_;
    foreach (@sitesets)  {
	ref($_) =~ /TFBS::SiteSet/ 
	    or $self->throw("Attempted to add an element ".
			    "that is not a TFBS::SiteSet object.");
	push @{$self->{site_list}},      @{ $_->{site_list} };
	push @{$self->{_iterator_list}}, @{ $_->{site_list} };
    }
    return $self;
}
	

    
sub sort  {
    my ($self) = @_;
    $self->{site_list} = [sort { $a->{start} <=> $b->{start} 
#				||
#				    $a->{end} <=> $b->{end}
				||
				$a->Matrix->{name} cmp $b->Matrix->{name}
			     }
			  @{$self->{site_list}}
			  ];
    $self->reset();
}
sub sort_reverse  {
    my ($self) = @_;
    $self->{site_list} = [sort { $b->{start} <=> $a->{start} 
#				||
#				    $a->{end} <=> $b->{end}
				||
				$b->Matrix->{name} cmp $a->Matrix->{name}
			     }
			  @{$self->{site_list}}
			  ];
    $self->reset();
}

sub sort_by_name  {
    my ($self) = @_;
    $self->{site_list} = [sort {$a->{name} cmp $b->{name} 
		   ||
		   $a->{start} <=> $b->{start} 
		   ||
		   $a->{end}   <=> $b->{end}	
	       }
	     @{$self->{site_list}}
	    ];
    $self->reset();
}

sub reset  {
    my ($self) = @_;
    @{$self->{_iterator_list}} = @{$self->{site_list}};
}


sub next_site  {
    my ($self) = @_;
    return shift @{$self->{_iterator_list}};
}

sub GFF  {
    use GFF::GeneFeatureSet;
    my ($self) = @_;
    $self->reset();
    my $GFFset = GFF::GeneFeatureSet->new(2);
    while (my $site = $self->next_site())  {
	$GFFset->addGeneFeature($site->GFF());
    }
    return $GFFset;
	
}



1;



require 5.005;
use ExtUtils::MakeMaker;


my $NAME     = 'TFBS';
my $DISTNAME = "TFBS";
my $VERSION  = "0.1";

get_sql_data();

WriteMakefile(
	      NAME	  => $NAME,
	      DISTNAME    => $DISTNAME,
	      VERSION     => $VERSION,
	      'dist'      => {  COMPRESS     => 'gzip -9f',
				SUFFIX       => '.gz', 
				DIST_DEFAULT => 'all tardist',
			    },
	      PREREQ_PM => { Bio::Root::RootI => 0.7,
		             PDL => 1.1,
		             GD  => 1.3,
			     File::Temp => 0.0,
			     IO::String => 0.0 
			    }, # e.g., Module::Name => 1.1
);

sub get_sql_data  {
    my $ans = "abc";
    do {
	print "Do you have write access to a MySQL database server? [n] ";
	$ans=<STDIN>; chomp $ans;
    } until $ans =~ /^y|n/i or $ans eq "";

    if (uc(substr($ans,0,1)) eq 'Y')  {
	print "\nOK, tell me more about it.\n\n";
	print "\tHost name : [localhost] ";	
	my $hostname = <STDIN>; chomp $hostname;
	$hostname = 'localhost' unless $hostname;

	print "\tUsername  : [none]     ";
	my $username = <STDIN>; chomp $username;
	$username = '' unless $username;
	
	print "\tPassword  : [none]     ";
	my $password = <STDIN>; chomp $password;
	$password = '' unless password;

	open FILE, ">t/MYSQLCONNECT" 
	    or die "Can't write to t/ directory, stopped";
	print FILE join("::", $hostname, $username, $password, " ");
	close FILE;
    }
    else {
	unlink "t/MYSQLCONNECT" if -e "t/MYSQLCONNECT";
    }
}

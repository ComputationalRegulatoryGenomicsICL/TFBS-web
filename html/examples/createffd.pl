use TFBS::DB::JASPAR2;
use TFBS::DB::FlatFileDir;
my $db = TFBS::DB::JASPAR2->connect("dbi:mysql:JASPAR2:forkhead","krivan",
				    "wk3003");
my $mset = $db->get_MatrixSet(-matrixtype=>"PFM");
my $mxit = $mset->Iterator;
my $i=0;
my $ffdb = TFBS::DB::FlatFileDir->create("examples/SAMPLE_FlatFileDir");
while (my $pfm = $mxit->next) {
    unless ($i++ % 7) {$ffdb->store_Matrix($pfm);}
}


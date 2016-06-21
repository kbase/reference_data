#!/usr/bin/env perl

use strict;
use FindBin qw($Bin);
use Getopt::Long::Descriptive;
use Data::Dumper;

my $solrServer = $ENV{KBASE_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";

my ($opt, $usage) = 
	describe_options(
		"%c %o",
		["division=s", "Division: bacteria | archaea | plant | fungi, multivalued, comma-seperated"],
		["source=s", "Source: genbank | refseq", {default => "refseq"}],
		["status=s", "Status: new | replaced | all", {default => "new"}],
		["format=s", "Format: gbf", { default => "gbf"}],
		["workspace=s", "Workspace", { default => "RefSeq_Genomes"}],
		["dir_path=s" , "ftp_path: path to the ftp directory for the reference genomes", { default => "genomes" }],
		[],
		["help|h", "Print usage message and exit"] );

print($usage->text), exit 0 if $opt->help;
print($usage->text), exit 1 unless $opt->division;
print($usage->text), exit 1 unless $opt->status;


my @divisions = split /,/, $opt->division;
my $source= $opt->source;
my $status = $opt->status;
my @formats= split /,/, $opt->format;
my $workspace = $opt->workspace;

my $dir = $opt->dir_path;
`mkdir $dir` unless (-d $dir);

my $count = 0;

foreach my $division (@divisions){

	open TAB , ">$dir/genome_summary_$division.txt";

	my @assemblies = getAssemblySummary($source, $division);

	foreach my $entry (@assemblies){

		chomp $entry;

		if ($entry=~/^#/) { #header
			my $header = $entry;
			$header =~s/^#/status\tdivision\t/;
			print TAB "$header\n";
			next;
		}

		my @attribs = split /\t/, $entry;

		my $assembly = {};
		$assembly->{accession} = $attribs[0];
		$assembly->{status} = $attribs[10];
		$assembly->{name} = $attribs[15];
		$assembly->{ftp_dir} = $attribs[19];
		$assembly->{file} = $assembly->{ftp_dir};
		$assembly->{file}=~s/.*\///; 
		
		($assembly->{id}, $assembly->{version}) = $assembly->{accession}=~/(.*)\.(\d+)$/;
		$assembly->{dir} = $assembly->{accession}."_".$assembly->{name};

		next if ($assembly->{status} eq "replaced"); 	

		print "Processing assembly: $assembly->{accession}\n";  
		
		my $assembly_status = checkAssemblyStatus($assembly);

		print TAB "$assembly_status\t$assembly->{genome_id}\t$division\t$entry\n";

		if ($assembly_status=~/(new|updated)/i){

			$count ++;

			my $genbankFile = getGenBankFile($assembly) if grep $_ eq "gbf", @formats;
			upload_genbank_file_to_ws_gto($workspace, $genbankFile);
			index_ws_gto_to_solr($workspace, $assembly->{accession});
		
		}else{

			# Current version already in KBase, check for annotation updates
			
		}

	last if $count>=5;

	}

}

close TAB;


sub getAssemblySummary {

	my ($source, $division) = @_;

	my $assembly_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/$source/$division/assembly_summary.txt";

	my @assemblies = `wget -q -O - $assembly_summary_url`;

	return @assemblies;

}


sub checkAssemblyStatus {

	my ($assembly) = @_;
	
	print "\tChecking status for assembly $assembly->{accession}: ";

	my $status;

  my $core = "/genomes";
  my $query = "/select?q=genome_id:".$assembly->{id}."*"; 
  my $fields = "&fl=genome_source,genome_id,genome_name";
  my $rows = "&rows=100";
  my $sort = "";
  my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	print "\n$solrQuery\n";

  my @records = `wget -q -O - "$solrQuery" | grep -v genome_name`;

	if (scalar @records == 0 ){
		$status = "New genome";
	}else{
		my ($genome_source, $genome_id, $genome_name) = split /\t/, @records[0];

		if ($genome_id eq $assembly->{accession}){
			$status = "Existing genome: current";
			$assembly->{genome_id} = $genome_id;
		}elsif ($genome_id =~/$assembly->{id}/){
			$status = "Existing genome: updated ";
			$assembly->{genome_id} = $genome_id;
		}else{
			$status = "Existing genome: status unknown";
			$assembly->{genome_id} = $genome_id;
		}
	}  

	print "$status\n";

	return $status;

}


sub getGenBankFile {

	my ($assembly) = @_;

	print "\tRetrieving GenBank file for $assembly->{accession}: ";

	my $url = "$assembly->{ftp_dir}/$assembly->{file}\_genomic.gbff.gz";
	my $outfile = "$dir/$assembly->{accession}.gbff.gz"; 
	
	my $err;	
	for (my $try=0; $try<5; $try++){
		$err = system("wget -q $url -O $outfile");
		last if $err == 0; 
	}

	my $md5 = `md5sum $outfile | perl -pe 's/ .*\n//'`;
	my $md5_expected = get_expected_md5($assembly);

	if ($err == 0 && $md5 eq $md5_expected){
			`gzip -df $outfile`;
			print "Success\n";
	}else{
			`rm $outfile`;
			print "Error => wget error code:$err\tMD5:$md5 vs $md5_expected\n";
	}

	$outfile=~s/\.gz$//;

	return $outfile;

}


sub get_expected_md5{

	my ($assembly) = @_;

	print "\tRetrieving md5checksums.txt file for $assembly->{accession}\n";

	my $url = "$assembly->{ftp_dir}/md5checksums.txt";
	
	my $row = `wget -q $url -O - | grep "genomic.gbff.gz"`;

	my $md5=$1 if $row=~/(\S*)\s.*genomic.gbff.gz/;

	return $md5;

}


sub upload_genbank_file_to_ws_gto {

	my ($workspace, $genbankFile) = @_;

	print "Uploading GenBank File $genbankFile to workspace as GTO\n";

	my $cmd = "time upload-genome -n --ci --workspace $workspace $genbankFile";

	print "$cmd\n";
	
	`$cmd`;		
	
}


sub index_ws_gto_to_solr {

	my ($workspace, $genome_id) = @_;

	print "Indexing genome $genome_id in Solr\n";

	my $cmd = "time $Bin/ws_genome_to_solr.pl --index yes --workspace $workspace --genome $genome_id";

	print "$cmd\n";

	`$cmd`;

}



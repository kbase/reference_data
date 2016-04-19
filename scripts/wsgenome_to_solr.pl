#!/usr/bin/env perl

use strict;
use FindBin qw($Bin);
use Getopt::Long::Descriptive;
use Data::Dumper;
use JSON;

my $json = JSON->new->allow_nonref;

# Define parameters 
my($opt, $usage) = describe_options("%c %o ",
			["genome=s" => 'WS genome object file'],
			["contigset=s" => 'WS contig set file'],
			["workspace=s" => "Workspace name", { default => 'RefSeq'}],
			["solr_url=s" => "Solr URL", { default => 'http://localhost:8005/solr'}],
			["index=s" => "Index genome in solr, yes|no", { default => 'yes'}],
			["help|h" => 'Show this help message'],
		);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->genome && $opt->contigset;

# Copy parameter values
my $ws_name = $opt->workspace;
my $ws_genome_name = $opt->genome;
my $ws_contigset_name = $opt->contigset;

my $genome_name = $ws_genome_name;
$genome_name =~s/(.gbff|.gbf|.gb)$//;

my $contigset_name = $ws_contigset_name;
$contigset_name =~s/(.gbff|.gbf|.gb)//;

#print "$genome_name\t$contigset_name\n";

#print parameters
#print "$ws_name:$genome_name:$contigset_name\n";

# Retrieve genome and contigset objects from the workspace

`ws-rename --workspace $ws_name $ws_genome_name $genome_name` unless $ws_genome_name eq $genome_name;
`ws-rename --workspace $ws_name $ws_contigset_name $contigset_name` unless $ws_contigset_name eq $contigset_name;

my $ws_genome_metadata  = `ws-get -w $ws_name $genome_name -m`;
my $ws_genome  = $json->decode(`ws-get -w $ws_name $genome_name`);
my $ws_contigset = $json->decode(`ws-get -w $ws_name $contigset_name`);

#print Dumper ($ws_genome_metadata, $ws_genome);

my @genome_metadata = split(/\n/, $ws_genome_metadata);
foreach my $metadata (@genome_metadata){
	my ($ws_genome_id) = $metadata=~/Object ID:(\d+)/ if $metadata=~/Object ID:(\d+)/;
}


=for
# Solr attributes for genome and genomic features
my @solr_keys = qw(object_id workspace_name object_type object_name genome_id feature_id genome_source genome_source_id feature_source_id protein_translation_length dna_sequence_length feature_type function gene_name aliases scientific_name scientific_name_sort genome_dna_size num_contigs num_cds complete domain taxonomy gc_content genome_publications feature_publications location_contig location_begin location_end location_strand locations roles subsystems subsystem_data protein_families annotations regulon_data atomic_regulons coexpressed_fids co_occurring_fids has_publications has_protein_families cs_db_version);
my @solr_genome_keys = qw(genome_id genome_source genome_source_id scientific_name scientific_name_sort genome_dna_size num_contigs num_cds complete domain taxonomy gc_content genome_publications has_publications);
my @solr_feature_keys = qw(feature_id feature_source_id protein_translation_length dna_sequence_length feature_type function gene_name aliases feature_publications location_contig location_begin location_end location_strand locations roles subsystems subsystem_data protein_families annotations regulon_data atomic_regulons coexpressed_fids co_occurring_fids has_protein_families);

#print Dumper  ([\@solr_keys, \@solr_genome_keys, \@solr_feature_keys]); 
=cut


my @solr_records;
my $record;
my $genome;
my @features;

# Prepare genome record for solr

$record->{workspace_name} = $ws_name; # KBasePublicRichGenomesV5
$record->{object_id} = $genome_name; #"kb|ws.".$ws_id.".obj."."$ws_genome_id"; # kb|ws.2869.obj.9837
$record->{object_name} = $genome_name; # kb|g.3397
$record->{object_type} = "KBaseGenomes.Genome-8.0"; # KBaseSearch.Genome-5.0 / KBaseGenomes.Genome-8.0

$record->{genome_id} = $genome_name; #$ws_genome->{id}; # kb|g.3397
$record->{genome_source} = "RefSeq"; # $ws_genome->{source}; # KBase Central Store
$record->{genome_source_id} = $ws_genome->{source_id}; # 83332.12

$record->{scientific_name} = $ws_genome->{scientific_name};
$record->{taxonomy} = $ws_genome->{taxonomy};
$record->{taxonomy} =~s/ *; */;/g;
#$genome->{tax_id} = $ws_genome->{tax_id};
$record->{domain} = $ws_genome->{domain};
$record->{gc_content} = $ws_genome->{gc_content};

$record->{genome_dna_size} = $ws_genome->{dna_size};
$record->{num_contigs} = $ws_genome->{num_contigs};
$record->{complete} = $ws_genome->{complete};

$record->{num_cds} = 0;
foreach my $feature (@{$ws_genome->{features}}){
	$record->{num_cds}++ if $feature->{type} = 'CDS'; 
}

#$genome->{genome_publications}=$ws_genome->{};
#$genome->{has_publications}=$ws_genome->{};

#print Dumper ($genome);
push @solr_records, $record;


# Prepare feature records for solr
foreach my $feature (@{$ws_genome->{features}}){

	my $record;

	$record->{workspace_name} = $ws_name; # KBasePublicRichGenomesV5
	$record->{object_id} = $feature->{id}; # kb|ws.2869.obj.9836/features/kb|g.3397.peg.3821
	$record->{object_name} = $feature->{id}; # kb|g.3397.featureset/features/kb|g.3397.peg.3821
	$record->{object_type} = "KBaseSearch.Feature"; # KBaseSearch.Feature

	$record->{genome_id} = $genome_name; #$ws_genome->{id}; # kb|g.3397
	$record->{genome_source} = "RefSeq"; # $ws_genome->{source}; # KBase Central Store
	$record->{genome_source_id} = $ws_genome->{source_id}; # 83332.12

	$record->{scientific_name} = $ws_genome->{scientific_name};
	$record->{taxonomy} = $ws_genome->{taxonomy};
	$record->{taxonomy} =~s/ *; */;/g;

	$record->{feature_type} = $feature->{type}; 
	$record->{feature_id} = $feature->{id}; 
	$record->{feature_source_id} = $feature->{id}; 
	$record->{gene_name} = $feature->{aliases}[0] unless $feature->{aliases}[0]=~/^(NP_|WP_|YP_|GI|GeneID)/i; # If not DBxref, must be gene name 
	$record->{function} = $feature->{function}; 
	$record->{aliases} = join(" :: ", @{$feature->{aliases}}); 
	
	$record->{location_contig} = $feature->{location}[0][0]; 
	$record->{location_begin} = $feature->{location}[0][1]; 
	$record->{location_end} = $feature->{location}[-1][1]+$feature->{location}[-1][3]; 
	$record->{location_strand} = $feature->{location}[0][2]; 
	$record->{locations} = $feature->{location}; 
	
	$record->{protein_translation_length} = $feature->{protein_translation_length}; 
	$record->{dna_sequence_length} = 0; 
	foreach my $location (@{$feature->{location}}){ $record->{dna_sequence_length}+= @{$location}[3]}; 

=for comment	
	$record->{roles} = $feature->{}; 
	$record->{subsystems} = $feature->{}; 
	$record->{subsystem_data} = $feature->{}; 
	$record->{protein_families} = $feature->{}; 
	$record->{annotations} = $feature->{annotations}; 
	$record->{regulon_data} = $feature->{}; 
	$record->{atomic_regulons} = $feature->{}; 
	$record->{coexpressed_fids} = $feature->{}; 
	$record->{co_occurring_fids} = $feature->{}; 
	$record->{has_protein_families} = $feature->{};	
	$record->{feature_publications} = $feature->{}; 
=cut

	push @solr_records, $record;

}


#print Dumper (\@solr_records);

my $genome_json = $json->pretty->encode(\@solr_records);
my $genome_file = "$genome_name.json";

open FH, ">$genome_file" or die "Cannot write to genome.json: $!";
print FH "$genome_json";
close FH;

`$Bin/post_solr_update.sh genomes $genome_file` if $opt->index=~/y|yes|true|1/i;


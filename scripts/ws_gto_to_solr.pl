#!/usr/bin/env perl

use strict;
use FindBin qw($Bin);
use Getopt::Long::Descriptive;
use Data::Dumper;
use JSON;

my $json = JSON->new->allow_nonref;

my @solr_records;

# Define parameters 
my($opt, $usage) = describe_options("%c %o ",
			["genome=s" => "Workspace genome object name", {default => 'GCF_000195955.2'}],
			["workspace=s" => "Workspace name", { default => 'test3'}],
			["solr_url=s" => "Solr URL", { default => 'http://localhost:8005/solr'}],
			["index=s" => "Index genome in solr, yes|no", { default => 'no'}],
			["help|h" => 'Show this help message'],
		);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->genome;

# Copy parameter values
my $ws_name = $opt->workspace;
my $ws_genome_name = $opt->genome;
my $genome_name = $ws_genome_name; ## tmp
my ($genome_source) = $ws_name=~/^(.*)_Genomes/;


# start building solr records

# Prepare genome record for solr

my $record;

my $ws_genome_metadata  = `ws-get -w $ws_name $ws_genome_name -m`;
my @genome_metadata = split(/\n/, $ws_genome_metadata);
foreach my $metadata (@genome_metadata){
	my ($ws_genome_id) = $metadata=~/Object ID:(\d+)/ if $metadata=~/Object ID:(\d+)/;
}
$record->{workspace_name} = $ws_name; # KBasePublicRichGenomesV5
$record->{object_id} = $ws_genome_name; #"kb|ws.".$ws_id.".obj."."$ws_genome_id"; # kb|ws.2869.obj.9837
$record->{object_name} = $ws_genome_name; # kb|g.3397
$record->{object_type} = "KBaseGenomes.Genome-8.0"; # KBaseSearch.Genome-5.0 / KBaseGenomes.Genome-8.0

# Get genome info
my $ws_genome  = $json->decode(`ws-get -w $ws_name $ws_genome_name`);
$record->{genome_id} = $ws_genome_name; #$ws_genome->{id}; # kb|g.3397
$record->{genome_source} = $genome_source; # $ws_genome->{external_source}; # KBase Central Store
$record->{genome_source_id} = $ws_genome->{external_source_id}; # 83332.12
$record->{num_cds} = $ws_genome->{counts_map}->{CDS};

# Get assembly info
my $ws_assembly = $json->decode(`ws-get $ws_genome->{assembly_ref}`);
$record->{genome_dna_size} = $ws_assembly->{dna_size};
$record->{num_contigs} = $ws_assembly->{num_contigs};
$record->{complete} = ""; #$ws_genome->{complete}; # type?? 
$record->{gc_content} = $ws_genome->{gc_content};

# Get taxon info
my $ws_taxon = $json->decode(`ws-get $ws_genome->{taxon_ref}`);
$record->{scientific_name} = $ws_taxon->{scientific_name};
$record->{taxonomy} = $ws_taxon->{scientific_lineage};
$record->{taxonomy} =~s/ *; */;/g;
#$record->{tax_id} = $ws_taxon->{taxonomy_id};
$record->{domain} = $ws_taxon->{domain};

#$genome->{genome_publications}=$ws_genome->{};
#$genome->{has_publications}=$ws_genome->{};

push @solr_records, $record;

#print Dumper(\@solr_records);


# Prepare feature records for solr

foreach my $feature_type (keys %{$ws_genome->{feature_container_references}}){

 my $container_ref = $ws_genome->{feature_container_references}->{$feature_type};
 my $ws_features = $json->decode(`ws-get $container_ref`);

 #print Dumper ($ws_features);

 foreach my $key (keys %{$ws_features->{features}}){

	my $feature = $ws_features->{features}->{$key};

	my $record;

	$record->{workspace_name} = $ws_name; # KBasePublicRichGenomesV5
	$record->{object_id} = $feature->{feature_id}; # kb|ws.2869.obj.9836/features/kb|g.3397.peg.3821
	$record->{object_name} = $feature->{feature_id}; # kb|g.3397.featureset/features/kb|g.3397.peg.3821
	$record->{object_type} = "KBaseSearch.Feature"; # KBaseSearch.Feature

	$record->{genome_id} = $ws_genome_name; #$ws_genome->{id}; # kb|g.3397
	$record->{genome_source} = $genome_source; # $ws_genome->{external_source}; # KBase Central Store
	$record->{genome_source_id} = $ws_genome->{external_source_id}; # 83332.12

	$record->{scientific_name} = $ws_taxon->{scientific_name};
	$record->{taxonomy} = $ws_taxon->{scientific_lineage};
	$record->{taxonomy} =~s/ *; */;/g;
	#$record->{tax_id} = $ws_taxon->{taxonomy_id};
	$record->{domain} = $ws_taxon->{domain};

	$record->{feature_type} = $feature->{type}; 
	$record->{feature_id} = $feature->{feature_id}; 
	$record->{feature_source_id} = $feature->{feature_id}; 
	$record->{function} = $feature->{function}; 
	
	# aliases
	my @aliases;
	foreach my $key (keys %{$feature->{aliases}}){
		push @aliases, "$feature->{aliases}->{$key}[0]:$key";
		$record->{gene_name} = $key if $feature->{aliases}->{$key}[0]=~/Genbank Gene/i;	 
	}
	$record->{aliases} = join(" :: ", @aliases);
	
	my $last_location = scalar @{$feature->{locations}};
	$record->{location_contig} = $feature->{locations}[0][0]; 
	$record->{location_begin} = $feature->{locations}[0][1]; 
	$record->{location_end} = $feature->{locations}[$last_location][1]+$feature->{location}[$last_location][3]; 
	$record->{location_strand} = $feature->{locations}[0][2]; 
	$record->{locations} = $json->pretty->encode($feature->{locations});
	$record->{locations} =~s/\s*//g; 
	$record->{locations} =~s/,\[\]//g; 

	$record->{protein_translation_length} = $feature->{protein_translation_length}; 
	$record->{dna_sequence_length} = $feature->{dna_sequence_length}; 

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

}


#print Dumper (\@solr_records);

my $genome_json = $json->pretty->encode(\@solr_records);
my $genome_file = "$genome_name.json";

open FH, ">$genome_file" or die "Cannot write to genome.json: $!";
print FH "$genome_json";
close FH;

`$Bin/post_solr_update.sh genomes $genome_file` if $opt->index=~/y|yes|true|1/i;


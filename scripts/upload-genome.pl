use Data::Dumper;
use strict;
use File::Basename;
use HTTP::Request::Common;
use Bio::KBase::userandjobstate::Client;
use Bio::KBase::Transform::Client;
use JSON::XS;
use LWP::UserAgent;
use File::Basename;
use Getopt::Long::Descriptive;
use Time::HiRes 'gettimeofday';
use Proc::ParallelLoop;
use Bio::KBase::AuthToken;
use IO::Handle;

$| = 1;

my($opt, $usage) = describe_options("%c %o data-file [data-file ...]",
				    ["timeout=s" => 'Job timeout', { default => 120 } ],
				    ["log-dir=s" => 'Logging directory'],
				    ["use-new-genome-type|n" => "Use new genome annotation type for uploaded data"],
				    ["ci" => "Use the CI infrastructure"],
				    ["base-url=s" => "Use the given base url" ],
				    ["parallel=i" => "Parallel threads", { default => 1 }],
				    ["shock_service_url=s" => "Shock URL"],
				    ["handle_service_url=s" => "Handle service url"],
				    ["ujs_service_url=s" => "UJS url"],
				    ["transform_service_url=s" => "Transform service url"],
				    ["workspace=s" => "Workspace name"],
				    ["help|h" => 'Show this help message'],
				    );

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV < 1;

my $target_type = 'KBaseGenomes.Genome';
if ($opt->use_new_genome_type)
{
    $target_type = "KBaseGenomeAnnotations.GenomeAnnotation";
}

my @data_files = @ARGV;

# --shock_service_url https://ci.kbase.us/services/shock-api/ --handle_service_url https://ci.kbase.us/services/handle_service/ --ujs_service_url https://ci.kbase.us/services/userandjobstate/

my $tobj = Bio::KBase::AuthToken->new();
my $token = $tobj->token;

my($user) = $token =~ /un=([^|]+)/;

$ENV{KB_AUTH_TOKEN} = $token;

my $ua = LWP::UserAgent->new;
my $json = JSON::XS->new->pretty(1);

my $ws = $opt->workspace;

my $base_url = "https://kbase.us/services";
if ($opt->ci)
{
    $base_url = "https://ci.kbase.us/services";
}
if ($opt->base_url)
{
    $base_url = $opt->base_url;
}

my $shock_url = $opt->shock_service_url // "$base_url/shock-api/";
my $ujs_url = $opt->ujs_service_url // "$base_url/userandjobstate/";
my $transform_url = $opt->transform_service_url // "$base_url/transform/";


my $ujs = Bio::KBase::userandjobstate::Client->new($ujs_url);
my $transform  = Bio::KBase::Transform::Client->new($transform_url);

my @auth_header = ("Authorization"  => "OAuth $token");

my @work;

my $i = 0;
for my $file (@data_files)
{
    if (open(F, "<", $file))
    {
	my $fl = <F>;
	close(F);
	if ($fl !~ /LOCUS/)
	{
	    warn "Skipping $file: no LOCUS line at start\n";
	    next;
	}
    }
    else
    {
	warn "Skipping $file: cannot open: $!\n";
	next;
    }

    my $ws_name = basename($file);
    push(@work, [$file, $ws_name, $i++]);
}

pareach \@work, sub {
    my($work) = @_;
    my($f, $n, $i) = @$work;
    upload_file($f, $n, $i);
}, { Max_Workers => $opt->parallel };


sub upload_file
{
    my($path, $ws_name, $idx) = @_;

    my $log_fh;
    if ($opt->log_dir)
    {
	my $logfile = $opt->log_dir . "/$ws_name";

	# print "logfile=$logfile\n";
	# system("ls", "-l", $logfile);
	if (-s $logfile)
	{
	    print STDERR "Skipping already-processed $path\n";
	    return;
	}
	open($log_fh, ">",  $logfile) or die "Cannot open $logfile: $!";
	print STDERR "Logging $path to $logfile\n";
	$log_fh->autoflush(1);
    }

    my $t1 = gettimeofday;

    my $nres = $ua->post("$shock_url/node",
			 @auth_header);
    
    if (!$nres->is_success)
    {
	die "shock post failed: " . $nres->status_line . " " . $nres->content . "\n";
    }
    
    my $ndata = $json->decode($nres->content);
    
    if ($ndata->{status} ne 200)
    {
	die "Bad status from shock node request\n" . Dumper($ndata);
    }
    my $node_id = $ndata->{data}->{id};
    print $log_fh "Got node id $node_id\n" if $log_fh;
    print "Got node id $node_id\n";
    
    my $req = HTTP::Request::Common::POST("$shock_url/node/$node_id",
					  @auth_header,
					  Content_Type => 'multipart/form-data',
					  Content => [upload => [$path]]);
    $req->method('PUT');
    
    my $upres = $ua->request($req);
    if (!$upres->is_success)
    {
	print $log_fh "Upload put failed: " . $upres->message . " " . $upres->content . "\n" if $log_fh;
	die "Upload put failed: " . $upres->message . " " . $upres->content . "\n";
    }
    my $updata = $json->decode($upres->content);
    
    my $t2 = gettimeofday;
    printf $log_fh "$path\t$ws_name\tupload_time\t%f\n", $t2 - $t1 if $log_fh;
    printf STDERR "$path\t$ws_name\tupload_time\t%f\n", $t2 - $t1;
    
    my $params = {
	external_type => 'Genbank.Genome',
	kbase_type => $target_type,
	workspace_name => $ws,
	object_name => $ws_name,
	optional_arguments => {
	},
	url_mapping => {
	    'Genbank.Genome' => "$shock_url/node/$node_id",
	},
    };
    my $tres = $transform->upload($params);
    
    my($awe_id, $job_id) = @$tres;
    print $log_fh "Transforming: awe=$awe_id job=$job_id\n" if $log_fh;
    print "Transforming: awe=$awe_id job=$job_id\n";

    my $t3 = gettimeofday;
    printf $log_fh "$path\t$ws_name\ttransform_start\t%f\n", $t3 - $t2 if $log_fh;
    printf STDERR "$path\t$ws_name\ttransform_start\t%f\n", $t3 - $t2;

    my $last_status = "awaiting_start";
    my $last_t = $t3;

    #
    # Wait for two minutes.
    #
    my $drop_dead;
    my $ok;

    while (!defined($drop_dead) || (time < $drop_dead))
    {
	my @res = $ujs->get_job_status($job_id);
#	print Dumper(\@res);
	
	my($last_update, $stage, $status, $progress, $est_complete, $complete, $error) = @res;

	if ($status =~ /Initializ/ && !defined($drop_dead))
	{
	    $drop_dead = time + $opt->timeout;
	    print "Set dropdead with status=$status\n";
	    print $log_fh "Set dropdead with status=$status\n" if $log_fh;
	}

	if ($status ne $last_status)
	{
	    my $t = gettimeofday;
	    printf $log_fh "$path\t$ws_name\tstatus=$last_status\t%f\n", $t - $last_t if $log_fh;
	    $path =~ s/\r/_r/g;
	    $ws_name =~ s/\r/_r/g;
	    $last_status =~ s/\r/_r/g;
	    printf "$path\t$ws_name\tstatus=$last_status\t%f\n", $t - $last_t;
	    $last_t = $t;
	    $last_status = $status;
	}
	
	if ($complete)
	{
	    print $log_fh "Complete: status=$status\n" if $log_fh;
	    print "Complete: status=$status\n";
	    $ok++;
	    last;
	}
	sleep 1;
    }

    my $det_err = $ujs->get_detailed_error($job_id);
    $det_err =~ s/\r//g;
    if ($det_err)
    {
	print $log_fh Dumper($det_err) if $log_fh;
	print Dumper($det_err);
    }

    if (!$ok)
    {
	print $log_fh "Job timed out\n" if $log_fh;
	print "Job timed out\n";
    }
    
    my $t4 = gettimeofday;
    printf $log_fh "$path\t$ws_name\ttransform\t%f\n", $t4 - $t3 if $log_fh;
    printf $log_fh "$path\t$ws_name\ttotal\t%f\n", $t4 - $t1 if $log_fh;
    printf STDERR "$path\t$ws_name\ttransform\t%f\n", $t4 - $t3;
    printf STDERR "$path\t$ws_name\ttotal\t%f\n", $t4 - $t1;
}

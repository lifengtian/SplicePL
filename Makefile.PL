# $Id$
use ExtUtils::MakeMaker;

eval "use Test::Manifest 1.14";

my $script_name = "SplicePL.pm";

WriteMakefile(
		'NAME'      => $script_name,
        'VERSION'   => '0.10',

		'EXE_FILES' =>  [ $script_name ],
            
		'PREREQ_PM' => {
                },

		'MAN1PODS'  => {
			$script_name => "\$(INST_MAN1DIR)/$script_name.1",
			},

		clean => { FILES => "*.bak $script_name-*" },
		);
        
1;

**coveragekit**
---------------

----------

coveragekit is a python package and associated number of scripts that can be used to calculate QC metrics needed for NGS sequencing (targeted panel, exome, or whole-genome) in a clinical context. There are a lot of coverage tools out there but they all seem to fall down somewhere in the QC process. 

In particular coveragekit aims to provide the following functionality:

 - (Relatively) fast bam parsing.
 - Structured JSON output.
 - Correct handling of irksome edge cases such as duplicates, overlapping mate pairs, insertion/deletions, and read clipping (more detail below).
 - Creation of indexed SQLite database with coverage metrics per gene to allow for fast subsequent exploration of coverage metrics for any set of genes.

**Installation**
----------------

coveragekit has been tested using Python 2.7.6+ and has the following dependencies:

 - pysam(0.8.4)

coveragekit can be installed like so:

    git clone git@gitlab.cs.washington.edu:cjhale/coveragekit.git
    cd coveragekit
    python setup.py install

**Usage**
---------

The easiest way to use coveragekit is by using the included coveragekit.py script from the command line. Currently it provides two functionalities:

 - bam - a way to parse bam files
 - db - an easy way to perform queries on a coveragekit SQLite database


bam
---

Running "coveragekit.py bam -h" will give the following result:

    Usage: coveragekit bam --bam sample.bam

    Options:
      -h, --help            show this help message and exit
      -b BAM, --bam=BAM     Input bam.
      -r REGIONS, --regions=REGIONS
                            Region file in bed format prepended with colon-
                            delimited descriptor ( eg 'reference:file.bed' ).
      -d DATABASES, --databases=DATABASES
                            Database files to build prepended with colon-delimited
                            descriptor to match region file ( eg
                            'reference:file.db' ).
      -w WINDOWSIZE, --windowSize=WINDOWSIZE
                            Processing window size [1000000].
      -t THREADS, --threads=THREADS
                            Number of processing threads.
      -l LEVELS, --levels=LEVELS
                            Comma-separated coverage levels for reporting
                            ['5,10,20,50,100'].
      --mq=MAPQ             Mapping quality cutoff [1].
      --genome              Calculate coverage for a genome [False].
      --allowdups           Count duplicate reads [False].
      --json=JSON           Output file for json doc.
      --txt=TXT             Output file for txt report.

 Some of these options are fairly self explanatory, some are less so. The easiest way to explain all the options is to present an example. If we wanted to assay coverage of an exome experiment we might run something like the following:

    python coveragekit.py bam \
      --bam exome.bam \
      --regions exome_target:exome_target.bed \
      --databases exome_target:exome_target_coverage.db \
      --windowSize 1000000 \
      --threads 4 \
      --levels 4,8,16 \
      --mq 20 \
      --json exome_target_coverage_results.json \
      --txt exome_target_coverage_results.txt ;

Such a call would cause coveragekit to read the exome.bam input file (make sure it's indexed!) and use the exome_target.bed file to calculate capture target coverage stats. You can specify as many bed files as you want through repeated invocation of the "--regions" flag, but each target must have a unique colon-delimited descriptor pre-pended to the bed file path. These descriptors are used in subsequent reporting so it helps if they are meaningful and consistent across exomes.

The "--databases" flag provides coveragekit with the output file path for an SQLite databases it generates. Like the "regions" argument, the "databases" argument can be specified multiple times. The pre-pended descriptor for database files must match one of "regions" file descriptors as there is a 1:1 relationship to a SQLite database and an input bed file. You do not have to create any databases, but if you do, there must be a paired region file.

The "windowSize" and "threads" arguments help tune performance. More threads are better, and the window size (which correlates to the amount of a bam file read at a time) does not matter unless you have very uneven distribution of target regions in the genome.

In the case above coveragekit will calculate coverage stats using the cutoffs specified by the "--levels" option (i.e. % covered at 4X, % covered at 8X, percent covered at 16X). These cutoffs define the database structure, so subsequent queries of the SQLite database will be limited to those coarse groupings.

The "--mq 20" argument above means only reads with a mapping quality >= 20 will be considered. If "--allowdups" is specified then duplicate reads will be counted (don't do this).

The "--json" and "--txt" files allow you to specify the paths for output in either json or tsv format. The details of these formats are below.

Finally, if you are processing a whole genome, you will want to specify "--genome" to force coveragekit to assay the depth of coverage at every basepair, rather than jumping from target to target. This mode is much slower than the default.

**Inputs**

It should be noted that you can pass in input bed files with an arbitrary columns, but coveragekit will only use the first three, plus the forth as a region name if it is present. Thus you should expect to input something like this:

    1	100	300	GeneA
    1	350	400	GeneA
    1	500	600	GeneB
    1	700	900	GeneC
    1	1100	1200	GeneC

Given this input coveragekit will automatically aggregate results for the first two rows into a record for "GeneA". "GeneB" will be assumed to be a single-region (exon) gene, and "GeneC" will be assigned the last two regions. If you do not give a gene name each region will be assigned an arbitrary numeric identifier. If you want regions to be grouped together, the must have identical identifiers in the fourth column. Also note that regions will not be grouped between multiple bed files.

**Outputs**

After running the bam command via coveragekit.py you should see a JSON output that looks something like this:

    {
        "allReads": 125062442,
        "inputBam": "exome.bam",
        "insertMean": 290.8247085006931,
        "insertSD": 90.74794926422638,
        "onTarget": {
            "exome_target": 73549166.0
        },
        "readsCounted": 121154783,
        "readsNotCounted": {
            "duplicate": 3839435,
            "mapquality": 0,
            "unmapped": 68224
        },
        "regionStats": {
            "exome_target": {
                "avgCoverage": 164.93806601866697,
                "coverageLevels": {
                    "0": 0.999194623137355,
                    "20": 0.986957068395623,
                    "100": 0.6817241269917957
                },
                "file": "exome_target.bed",
                "length": 34359070,
                "numRegions": 19241
            }
        },
        "version": "1.1.0"
    }
      
In this output there are details about the bam file that was parsed and two nested data structures pertaining to any bed regions specified on the command line, the "onTarget" and "regionStats" fields. For each of these members, the statistics pertaining to a particular bed file are keyed by the descriptor passed on the command line. All of the values in the root level of the output relate to read counts with the exception of the "insertMean" and "insertSD" fields which correspond to the length of the sequencing library inserts in bp. In the "regionStats" object, the numbers should be self explanatory with the exception of those within the "coverageLevels" member, which correspond to the ratio of base pairs within a given target covered at a level that corresponds with those levels passed at the command line.

The output specified by the "--txt" flag is a simple text formatted document that essential mimics the JSON output while being slightly more human readable.


db
--

Running "python coveragekit.py db -h" brings up the following:

    Usage: coveragekit db --db sample.db [ options ]
    
    Options:
      -h, --help            show this help message and exit
      -d DB, --db=DB        Input database.
      --geneList=GENELIST   Comma-separated gene list.
      --geneListFile=GENELISTFILE
                            File with newline-separated gene list.
      --levelsMin=LEVELSMIN
                            Comma-separated list of minimum percents at X coverage
                            with colon delimited coverage level prepended ( eg
                            '5:99,10:95,20:90' ).
      --levelsMax=LEVELSMAX
                            Comma-separated list of maximum percents at X coverage
                            with colon delimited coverage level prepended ( eg
                            '5:99,10:95,20:90' ).
      --coverageMin=COVERAGEMIN
                            Minimum average coverage.
      --coverageMax=COVERAGEMAX
                            Maximum average coverage.
      --reportRegions       Report regions with coverage of intersest as JSON
                            stings.
      --json=JSON           Output JSON file.
      --tsv=TSV, --txt=TSV  Output tsv file.


Here the important thing to remember is that the regions in the database you are going to query were defined by the input bed given to the "coveragekit.py bam" call and the coverage levels you can assay were defined by the "--levels" of that bam call. Therefore when using the "--levelsMin" or "--levelsMax" options the levels pre-pended to each percent must have been specified in the bam call, and any genes specified with the "--geneList" or "--geneListFile" options had to have been present and named consistently in the capture target bed files used as input.

Additionally, it is important to point out that the "levelsMin" option returns all regions with >= the specified percentage coverage at a given level whereas "levelMax" returns all regions with < the specified percentage coverage at a given level.

Thus, if you wanted to find if any of the genes in a list of 5 genes were covered at least 90% at 4X and less than 100% at 8X you would use the following:

    python coveragekit.py db \
      --db exome_target_coverage.db \
      --geneList GeneA,GeneB,GeneC,GeneD,GeneE \
      --levelsMin 4:90 \
      --levelsMax 8:100 \
      --reportRegions \
      --json exome_target_coverage.query_result.json \
      --txt exome_target_coverage.query_result.txt ;

Using the command above you may see the following result in a JSON file:

    {
        "meta": {
            "coverageSource": "exome.bam",
            "dbLevels": [
                0,
                4,
                8,
    			16
            ],
            "dbSource": "exome_target_coverage.db",
            "queryResultNum": 2,
            "queryString": "SELECT * FROM regions WHERE id IN (\"GeneA\",\"GeneB\",\"GeneC\",\"GeneD\",\"GeneE\") AND percent4X >= 0.9 AND percent8X < 1.0",
            "regionSource": "exome_capture.bed",
            "version": "1.1.0"
        },
        "queryResults": [
            {
                "coverage": 16.54343654,
                "coverageRegions": {
                    "greaterOrEqual": {
    					"4": [
                            "12:9247568-9247581",
                            "12:9262909-9262930"
                        ]
    				},
                    "lessThan": {
                        "8": [
                            "12:9247568-9247581",
                            "12:9262909-9262930"
                        ]
                    }
                },
                "id": "GeneB",
                "percentGreaterOrEqual": {
                    "0": 1.0,
                    "4": 0.9923163841807909,
                    "8": 0.4894915254237288,
    				"16": 0.1004850980234907
                },
                "position": "12:9247568-9262930"
            },
            {
                "coverage": 22.42236423,
                "coverageRegions": {
                    "greaterOrEqual": {
                        "4": [
                            "4:170991738-170991740"
                        ]
                    },
                    "lessThan": {
                        "8": [
                            "4:170991738-170991740"
                        ]
                    }
                },
                "id": "GeneE",
                "percentGreaterOrEqual": {
                    "0": 1.0,
                    "4": 0.9813953488372092,
                    "8": 0.4961240310077519,
    				"16": 0.0780934098098675
                },
                "position": "4:170991738-170991740"
            }
    	]
    }

In the JSON output there are two main branches, the "meta" section and the "queryResults" section. The "meta" section gives you information about the database queried, version of coveragekit used, and the actual SQL syntax used in the query.

The "queryResults" section is a list of genes (or numeric regions if gene names were not specified) returned by the query. Most of the statistics for a gene are aggregated across the entire gene, with the exception of the "coverageRegions" section which is only returned if "reportRegions" is specified in the command. In this section, the exact regions of a gene that fell above or below specified cutoffs are given. These can be useful in pinpointing exactly which bases failed to cross a coverage threshold.

The tsv output is essentially a representation of the JSON output with each row representing a gene or region.

> Written with [StackEdit](https://stackedit.io/).


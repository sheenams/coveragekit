#!/usr/bin/env python

import sys, optparse, json, datetime, logging, math
import pysam
import coveragekit.utils.region as covregion
import coveragekit.utils.db as covdb
import coveragekit.utils.bed as covbed
from coveragekit.utils.bam import BamReader,BamReaderAggregate,ProcessingRegionGenerator

from multiprocessing import Pool

from coveragekit.version import __version__

def _readBamRegion(inputs):
    bam = inputs[0]
    levels = inputs[1]
    regions = inputs[2]
    mapq = inputs[3]
    dups = inputs[4]
    genome = inputs[5]
    startTime = inputs[6]
    totalRegions = inputs[7]
    
    logger = logging.getLogger("bam reader thread")
    logger.setLevel(logging.INFO)
    
    updateSize = math.floor(totalRegions/float(100))
    if (regions[0].index % updateSize == 0) and (regions[0].index > 0):
        now = datetime.datetime.now()
        timeDiff = now - startTime
        remainingRegions = totalRegions - (regions[0].index +1)
        logger.info("Processing region {}\tTime elapsed - {:.2f}m\tTime remaining - {:.2f}m".format(regions[0].index, (timeDiff.seconds/60.0), ((timeDiff.seconds/float(regions[0].index+1))*remainingRegions)/60.0))
    bamRegion = BamReader(bam, regions, levels, mapq, dups, genome)
    bamRegion.read()
    return bamRegion.report()

def bam(bamInput, regions, databases, levels, windowSize, threads, mapq, dups, genome):
    '''Returns a dict containing coverage data information for a given bam file.
    
    :param bamInput: file path for bam file
    :type bamInput: str
    :param regions: Dict of regions to assay coverage over with key:value pairs of region descriptor:region file path
    :type regions: dict
    :param levels: List of integers corresponding to levels of coverage to consider
    :type level: list
    :param windowSize: Size of bam chunk to be considered by a bam reader
    :type windowSize: int
    :param threads: Number of BamReader processes to launch
    :type threads: int
    :param mapq: Minimum mapping quality score to make a read eligible for 
    :type mapq: int
    :param dups: Boolean indicating whether duplicate reads should be considered 
    :type dups: bool
    :param genome: Boolean indicating whether bam file should have genome-level coverage considered 
    :type genome: bool
    
    :rtype: dict
    
    '''
    
    # Set up logging
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger("coveragekit bam")
    logger.setLevel(logging.INFO)
    
    startTime = datetime.datetime.now()
    
    # Get a list of regionSets
    regionSets = regions.keys()
    logger.info("Preparing to read from {} input region files".format(len(regionSets)))
    
    # Parse region files
    processingRegionGenerator = ProcessingRegionGenerator(bamInput,windowSize)
    
    if len(regions) > 0:
        for descriptor,bedFile in regions.items():
            for bedRegion in covbed.bedToRegions(descriptor,bedFile):
                processingRegionGenerator.addRegion(bedRegion)
    
    # Create region set objects to aggregate all the stats over given capture or gene sets
    logger.info("Creating processing regions using specified window size and input regions".format(len(regionSets)))
    regionSetAggregators = {}
    for descriptor in regionSets:
        regionSetAggregators[descriptor] = covregion.RegionSet(descriptor, levels)
    
    # Get bam header
    #bamFile = pysam.AlignmentFile(bamInput, 'rb')
    #processingRegions = _defineProcessingRegions(bamFile.header, parsedRegions, windowSize)
    #bamFile.close()
    

    # Launch bam reading threads
    bamJobs = []
    for r in processingRegionGenerator.returnProcessingRegion():
        bamJobs.append((bamInput, tuple(levels), r, mapq, dups, genome, startTime))
    
    # The following loop is simply for run-time estimation in the multiprocessing forks
    for i in range(0,len(bamJobs)):
        bamJobs[i] += (len(bamJobs),)    
    
    logger.info("Total regions to process: {}".format(len(bamJobs)))
    
    bamWorkers = Pool(processes = threads)
    results = bamWorkers.map(_readBamRegion, bamJobs, chunksize=1)
    
    # Uncomment the follow for debugging purposes
    #results = []
    #for curJob in bamJobs[321:]:
    #    results.append(_readBamRegion(curJob))
    #results = [_readBamRegion(bamJobs[0])]
    
    # Now we parse the results for each chunk of alignment data
    bamAggregator = BamReaderAggregate(regionSets)
    for chunk in results:
        # Aggregate stats for the bam in question
        bamAggregator.add(chunk)
        
        # Add subregions to region aggregator objects
        for subRegionResult in chunk[8]:
            regionSetAggregators[subRegionResult[0].regionSet].add(subRegionResult[0],subRegionResult[1])
    
    # Reporting time
    report = bamAggregator.report(bamInput, genome)
    
    # The following supplements the BamReaderAggregate report with region reports
    report["regionStats"] = {}
    for descriptor in regionSets:
        regionReport = regionSetAggregators[descriptor].report()
        regionSetName = regionReport["name"]
        
        del regionReport["name"]
        report["regionStats"][regionSetName] = regionReport
        report["regionStats"][regionSetName]["file"] = regions[descriptor]
    
    # Create coverage databases
    for databaseKey,databaseFile in databases.items():
        referenceDB = covdb.CoverageDB(databaseFile,
                                       regionsource = regions[databaseKey],
                                       coveragesource = bamInput,
                                       levels = regionSetAggregators[databaseKey].levels,
                                       mapq = mapq,
                                       dups = dups,
                                       totalCoverage = bamAggregator.totalCoverage,
                                       overwrite = True)
        referenceDB.insertRegionSet(regionSetAggregators[databaseKey])
    
    logger.info("Finished.")
    
    return report
    
def report(data, jsonOut = None, txtOut = None):
    if txtOut:
        with open(txtOut, "w") as txtFH:
            txtFH.write("coveragekit bam (v{}) -- text report".format(data["version"]))
            txtFH.write("\n\nInput BAM file:\t{}\n".format(data["inputBam"]))
            txtFH.write("Text report file:\t{}\n".format(txtOut))
            if jsonOut:
                txtFH.write("JSON report file:\t{}\n".format(jsonOut))
            txtFH.write("\nTotal reads:\t{}\n".format(data["allReads"]))   
            txtFH.write("Number of reads counted:\t{}\n".format(data["readsCounted"]))
            txtFH.write("Number of reads not counted:\n")
            for key,value in data["readsNotCounted"].items():
                txtFH.write("\t{}:\t{:3.2f}%\t({})\n".format(key,(value/float(data["allReads"]))*100,value))
            txtFH.write("Average insert size estimate:\t{}\n".format(data["insertMean"]))
            txtFH.write("Insert size standard deviation estimate:\t{}\n".format(data["insertSD"]))
            
            if "genome" in data.keys():
                txtFH.write("Average genome-wide coverage:\t{}\n".format(data["genome"]["avgCoverage"]))
            
            txtFH.write("On target percentages:\n")
            for key,value in data["onTarget"].items():
                txtFH.write("\t{}:\t{:3.2f}%\n".format(key,((value/float(data["readsCounted"]))*100)))
            txtFH.write("Region stats:\n")
            for regionNames,stats in data["regionStats"].items():
                txtFH.write("\t{}:\n".format(regionNames))
                txtFH.write("\t\tRegion file:\t{}\n".format(stats["file"]))
                txtFH.write("\t\tNumber of regions:\t{}\n".format(stats["numRegions"]))
                txtFH.write("\t\tLength:\t{}\n".format(stats["length"]))
                txtFH.write("\t\tAverage Coverage:\t{}\n".format(stats["avgCoverage"]))
                txtFH.write("\t\tPercent at X coverage or greater:\n")
                for key,value in stats["coverageLevels"].items():
                    txtFH.write("\t\t\t{}X:\t{:3.2f}\n".format(key,(value*100)))
        
    if jsonOut:
        with open(jsonOut, "w") as jsonFH:
            jsonFH.write(json.dumps(data, indent=4, sort_keys=True))

def run(inputArgs):
    usage = "%prog --bam sample.bam"
    parser = optparse.OptionParser(usage=usage, prog = "coveragekit bam")
    parser.add_option("-b","--bam", dest="bam", help="Input bam.", default="")
    parser.add_option("-r","--regions", action="append", dest="regions", help="Region file in bed format prepended with colon-delimited descriptor ( eg 'reference:file.bed' ).", default=[])
    parser.add_option("-d","--databases", action="append", dest="databases", help="Database files to build prepended with colon-delimited descriptor to match region file ( eg 'reference:file.db' ).", default=[])
    parser.add_option("-w","--windowSize", type="int", dest="windowSize", help="Processing window size [1000000].", default=1000000)
    parser.add_option("-t","--threads", type="int", dest="threads", help="Number of processing threads.", default=1)
    parser.add_option("-l","--levels", type="string", dest="levels", help="Comma-separated coverage levels for reporting ['5,10,20,50,100'].", default="5,10,20,50,100")
    parser.add_option("--mq", type="int", dest="mapq", help="Mapping quality cutoff [1].", default=1)
    parser.add_option("--genome", action="store_true", dest="genome", help="Calculate coverage for a genome [False].", default=False)
    parser.add_option("--allowdups", action="store_true", dest="dups", help="Count duplicate reads [False].", default=False)
    parser.add_option("--json", type="string", dest="json", help="Output file for json doc.", default=None)
    parser.add_option("--txt", type="string", dest="txt", help="Output file for txt report.", default=None)
    (options, args) = parser.parse_args(inputArgs)

    # Bam file is required as well as one output
    if len(options.bam) == 0: parser.error("Missing bam sample, use --bam or -b.")
    if (options.json is None) and (options.txt is None): parser.error("Must specify an output with --json or --txt")
    
    # Multiple region files can be submitted
    regions = {}
    if len(options.regions) > 0:
        for curRegion in options.regions:
            descriptorSplit = curRegion.split(":",1)
            if len(descriptorSplit) != 2:
                parser.error("Region files must have colon-delimited descriptor prepended.")
            regions[descriptorSplit[0]] = descriptorSplit[1]
    
    # Multiple database files can be submitted
    databases = {}
    if len(options.databases) > 0:
        for curDB in options.databases:
            descriptorSplit = curDB.split(":",1)
            if len(descriptorSplit) != 2:
                parser.error("Database files must have colon-delimited descriptor prepended.")
            databases[descriptorSplit[0]] = descriptorSplit[1]
        if len(set(databases.keys()).difference(set(regions.keys()))) > 0:
            parser.error("Database descriptors must match colon-delimited region file descriptors.")
    # Convert string of levels into sorted list of levels
    levels = []
    for i in options.levels.split(','):
        levels.append(int(i))
    levels.sort()
    coverageReport = bam(options.bam, regions, databases, levels, options.windowSize, options.threads, options.mapq, options.dups, options.genome)
    report(coverageReport, options.json, options.txt)


if __name__ == '__main__':
    run(sys.argv[1:])
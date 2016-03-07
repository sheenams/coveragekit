import pysam, logging, math
import coveragekit.utils.levels as levelkit
import coveragekit.utils.regioncaller as regioncaller
import coveragekit.utils.region as regionkit

from coveragekit.version import __version__

class BamRegion(object):
    ''' Class that extends the :class:`Region` class by adding callers to the :class:`CoverageLevel` class, and the onTarget attribute which keeps track of on-target reads.

    '''

    def add(self, pos, depth):
        '''Makes a call to the add method of the coverageLevel attribute, adding the depth at a given position to the coverage statistics for a region.
        
        :param pos: Chromosome coordinate (bp position).
        :type pos: int
        :param depth: Depth of non-redundant, quality coverage
        :type depth: int
        
        '''
        self.coverageLevel.add(pos, depth)
    
    def addOverlap(self, pos, read):
        '''Adds a read to the onTarget attribute.
        
        :param pos: Chromosome coordinate (bp position).
        :type pos: int
        :param read: Read id
        :type read: str
        
        '''
        if pos >= self.region.start: # This if statement ensures against edge-case where read alignment end doesn't cover a region due to clipping
            self.onTarget.update([read])

    def report(self, ):
        '''Returns a tuple with basic coverage metrics for a BamRegion.
        
        :returns: Tuple with format of (region object describing the bam region, set of read names overlapping with the region, CoverageLevel report)        
        :rtype: tuple
        
        '''
        return (self.region, self.onTarget, self.coverageLevel.report())

    def __init__(self, region, levels):
        '''Initializer for BamRegion class.
        
        :param region: The :class:`Region` used to define :class:`BamRegion`
        :type region: Region
        :param levels: List of levels to consider for metric levels
        :type regions: list
        
        '''
        self.region = region
        self.onTarget = set()
        self.coverageLevel = levelkit.CoverageLevel(self.region.start, self.region.stop, levels)

class BamReader(object):
    '''Class that reads part of a bam file and report backs coverage stats. The largest part that can be read is a chromosome / contig  listed in bam header.'''

    def read(self, ):
        '''Initiates a read of the bam file in the regions specified by the class attributes. This methods really consists of two sections.
        In the first part reads are parsed from the bam file and depending on user input (mapping quality cutoff, duplicates allowed),
        they are used to create a pileup and also used to calculate percent overlap with target files.
        
        In the second part, the pileup is walked over and used to calculate depth of coverage.
        
        '''
        coverage = [0 for x in range(self.region.length)]
        
        # Set cutoffs for updating region caller
        sortedSubregionStarts = list(reversed(sorted(self.subregionStarts.keys())))
        sortedSubregionStops = list(reversed(sorted(self.subregionStops.keys())))
        if len(sortedSubregionStarts) > 0:
            subRegionBasement = sortedSubregionStarts.pop()
        else:
            subRegionBasement = float("inf")
            
        if len(sortedSubregionStops) > 0:
            subRegionCeiling = sortedSubregionStops.pop()
        else:
            subRegionCeiling = float("inf")
        
        overlapRegionCaller = regioncaller.RegionCaller()
        overlapRegionCaller["_self"] = self.subregions[0].addOverlap
        
        # The readTracker will keep track of insert lengths
        readTracker = {}        
        
        # Iterate over bam reads
        chunkCount = 0
        for bamRead in self.bamReads:
            if ((not bamRead.is_duplicate) or (self.allowdups)) and (bamRead.mapping_quality >= self.qualityCutoff) and (not bamRead.is_unmapped):
                # Read names in bam format don't necessarily distinguish between 1st or second read in pair, so we make this explicit
                if bamRead.is_read1:
                    readName = "{}.1".format(bamRead.query_name)
                else:
                    readName = "{}.2".format(bamRead.query_name)
                if bamRead.reference_start < self.region.start:
                    readStart = self.region.start
                    
                    # Want to keep track of reads hanging off of this chunk so that we don't double count
                    self.firstColumn.append(readName)
                else:
                    readStart = bamRead.reference_start
                
                if bamRead.reference_end > self.region.stop:
                    readStop = self.region.stop
                    
                    # Want to keep track of reads hanging off of this chunk so that we don't double count
                    self.lastColumn.append(readName)
                else:
                    readStop = bamRead.reference_end
                 
                # Coverage assessment using cigar string to figure out covered regions
                coveragePos = readStart
                insertLength = 0
                cigarTuples = list(reversed(bamRead.cigartuples))
                while len(cigarTuples) > 0:
                    cigar = cigarTuples.pop()
                    if (cigar[0] in [0,7,8,2,3]): # Alignment match, sequence match, sequence mismatch - all of these add to length of insert as well coverage profile. Deletion or skip handled below
                        
                        # If the aligned portion of the read extends past the start of the paired alignment the read is overlapping and we don't count this towards coverage or insert length
                        if ((cigar[1] + coveragePos) >= bamRead.next_reference_start) and (bamRead.is_proper_pair) and (bamRead.template_length >= 0):
                            endPoint = bamRead.next_reference_start - coveragePos
                            cigarTuples = [] # Causes loop to exit after this operation
                        else:
                            endPoint = cigar[1]
                        
                        # Take care of situation where a read spans a chunk
                        if (coveragePos + endPoint) > self.region.stop:
                            endPoint = self.region.stop - coveragePos
                        
                        # Increment coverage array
                        for c in range((coveragePos - self.region.start),((coveragePos - self.region.start) + endPoint)):
                            coverage[c] += 1
                            
                        coveragePos += endPoint
                        
                        # Increase insert length unless there is a deletion from reference or skipped reference. Counts towards coverage profile but not insert length
                        if cigar[0] not in [2,3]:
                            insertLength += endPoint
                        
                    elif (cigar[0] == 1): # Insertion to reference. Does not count towards insert length, but not coverage profile. Maybe this should include soft clipping
                        insertLength += cigar[1]
                        
                    #elif (cigar[0] in [4,5]): # Soft or hard clipping - Neither count towards insert length or coverage profile - it could be argued that soft clipping should
                    #    pass
                        
                # Update the overlap event handler
                while readStop >= subRegionBasement:
                    for i in self.subregionStarts[subRegionBasement]:
                        overlapRegionCaller[self.subregions[i+1].region.index] = self.subregions[i+1].addOverlap
                    if len(sortedSubregionStarts) > 0:
                        subRegionBasement = sortedSubregionStarts.pop()
                    else:
                        subRegionBasement = float("inf")
                        break
                while readStart >= subRegionCeiling:
                    for i in self.subregionStops[subRegionCeiling]:
                        del overlapRegionCaller[self.subregions[i+1].region.index]
                    if len(sortedSubregionStops) > 0:
                        subRegionCeiling = sortedSubregionStops.pop()
                    else:
                        subRegionCeiling = float("inf")
                        break
                overlapRegionCaller(readStop-1,readName)
                
                # Calculate insert size
                if (bamRead.is_proper_pair):
                    if bamRead.query_name in readTracker:
                        insertLength += readTracker[bamRead.query_name]
                        del readTracker[bamRead.query_name]
                        #if insertLength > 10000:
                        #    print bamRead
                        self.insertLengths.append(insertLength)
                    else:
                        readTracker[bamRead.query_name] = insertLength + (bamRead.next_reference_start - coveragePos)
                chunkCount += 1
            elif bamRead.is_unmapped:
                self.uncountedMetrics["unmapped"] += 1
            elif bamRead.is_duplicate and not self.allowdups:
                self.uncountedMetrics["duplicate"] += 1
            elif bamRead.mapping_quality < self.qualityCutoff:
                self.uncountedMetrics["mapquality"] += 1
        self.logger.debug(chunkCount)
        
        if (len(self.subregions) > 1) or (self.genome == True):
            # Reset cutoffs for updating new region caller
            sortedSubregionStarts = list(reversed(sorted(self.subregionStarts.keys())))
            sortedSubregionStops = list(reversed(sorted(self.subregionStops.keys())))  
            if len(sortedSubregionStarts) > 0:
                subRegionBasement = sortedSubregionStarts.pop()
            else:
                subRegionBasement = float("inf")
                
            if len(sortedSubregionStops) > 0:
                subRegionCeiling = sortedSubregionStops.pop()
            else:
                subRegionCeiling = float("inf")
            
            coverageRegionCaller = regioncaller.RegionCaller()
            coverageRegionCaller["_self"] = self.subregions[0].add
            
            # Iterate over the coverage profile
            pos = 0
            while pos < len(coverage):
                depth = coverage[pos]
                # Update the event handler
                while (pos + self.region.start) >= subRegionBasement: # If the current position is after or equal to the first base of the next set of regions in the caller, see if we need to add more regions to the caller
                    for i in self.subregionStarts[subRegionBasement]:
                        coverageRegionCaller[self.subregions[i+1].region.index] = self.subregions[i+1].add
                    if len(sortedSubregionStarts) > 0:
                        subRegionBasement = sortedSubregionStarts.pop()
                    else:
                        subRegionBasement = float("inf")
                        break
                        
                while (pos + self.region.start) >= subRegionCeiling: # If the current position is after or equal to the last base of a set of regions in the caller, delete those regions from the caller
                    for i in self.subregionStops[subRegionCeiling]:
                        del coverageRegionCaller[self.subregions[i+1].region.index]
                    if len(sortedSubregionStops) > 0:
                        subRegionCeiling = sortedSubregionStops.pop()
                    else:
                        subRegionCeiling = float("inf")
                        break
                
                if (len(coverageRegionCaller) > 1) or (self.genome == True):
                    # Make a call to the event handler
                    coverageRegionCaller((pos + self.region.start), depth)
                    pos += 1
                elif subRegionBasement < float("inf"):
                    pos = subRegionBasement - self.region.start
                else:
                    break
        
        self.readFinished = True

    def report(self, ):
        '''Returns a tuple summarizing coverage statistics for the region of the bam file read by this reader in the following format:
        
            (:class:`Region` object for this chunk,
            int number of reads,
            dict of on target by region set,
            :class:`CoverageLevel` report for this chunk,
            dict of reads in first column,
            dict of reads in last column,
            dict of uncounted stats,
            list of insert sizes,
            [(:class:`Region` object for subregion1, :class:`BamRegion` report for subregion1),...])
        
        '''
        
        # First get the stats for the super region or chunk itself
        chunkTotal = self.subregions[0].report()
        
        # Create sets for first and last column read names as well as dicts that we will eventually return
        fSet = set(self.firstColumn)
        lSet = set(self.lastColumn)
        fDict = {key : [] for key in self.firstColumn}
        lDict = {key : [] for key in self.lastColumn}
        
        # Then get the stats for all of the sub-regions
        subRegionStats = []
        onTargetSets = {}
        onTarget = {}
        if len(self.subregions) > 1:
            for subregion in self.subregions[1:]:
                subregionReport = subregion.report()
                subregionInfo = subregionReport[0]                
                subRegionStats.append((subregionReport[0], subregionReport[2]))
                
                # Do on target aggregation by regionSet
                if subregionInfo.regionSet in onTargetSets:
                    onTargetSets[subregionInfo.regionSet].update(subregionReport[1])
                else:
                    onTargetSets[subregionInfo.regionSet] = set(subregionReport[1])
        
        # Get on-target numbers per region set
        for r,o in onTargetSets.items():
            onTarget[r] = len(o)
            
            # Append region to first and last column reads
            for n in fSet.intersection(o):
                fDict[n].append(r)
            for n in lSet.intersection(o):
                lDict[n].append(r)
            
        
        # Make final report tuple
        report = (chunkTotal[0], len(chunkTotal[1]), onTarget, chunkTotal[2], fDict, lDict, self.uncountedMetrics, self.insertLengths, subRegionStats)
        return report

    def __init__(self, bam, region, levels, qualityCutoff = 1, allowdups = False, genome = False):
        '''Returns a tuple summarizing coverage statistics for the region of the bam file read by this reader in the following format:
        
        (region object for this chunk,
        integer number of reads,
        dict of on target by region set,
        coverage level report for this chunk,
        dict of reads in first column,
        dict of reads in last column,
        dict of uncounted stats,
        list of insert sizes,
        [(subregion region object1, subregion coverage report1),...])
        
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
        
        self.logger = logging.getLogger("bam reader")
        self.logger.setLevel(logging.INFO)
        self.readFinished = False
        self.region = region[0]
        self.firstColumn = []
        self.lastColumn = []
        self.insertLengths = []
        self.qualityCutoff = qualityCutoff
        self.allowdups = allowdups
        self.genome = genome
        self.logger.debug(self.genome)
        self.uncountedMetrics = {"unmapped" : 0, "duplicate" : 0, "mapquality": 0}
        
        # Get reads from the current chunk
        bamfh = pysam.AlignmentFile(bam, 'rb')
        self.bamReads = bamfh.fetch(reference=self.region.chrom, start=self.region.start, end=self.region.stop)
        
        # Make sure the pileup uses the right chromosome nomenclature, and then strip out that stupid "chr" if it's in there
        if self.region.chrom.startswith("chr"):
            self.region.chrom = self.region.chrom[3:]
        
        self.subregions = []
        self.subregionStarts = {}
        self.subregionStops = {}
        
        self.regionCaller = regioncaller.RegionCaller()

        # Create subregions: subregion[0] is actually the super region or chunk being interrogated by this BamReader object
        self.subregions.append(BamRegion(self.region, levels))
        #self.regionCaller[-1] = self._smartAdderGen(self.subregions[-1].add, self.regionCaller, -1, self.region.stop)
        #self.regionCaller["_self"] = self.subregions[-1].add
        
        for i,r in enumerate(region[1]):
            self.subregions.append(BamRegion(r, levels))
            if r.start in self.subregionStarts:
                self.subregionStarts[r.start].append(i)
            else:
                self.subregionStarts[r.start] = [i]
            
            if r.stop in self.subregionStops:
                self.subregionStops[r.stop].append(i)
            else:
                self.subregionStops[r.stop] = [i]
                
class BamReaderAggregate(object):
    
    def add(self, results):
        
        # Collect the stats for the chunk as a whole
        resultsRegion = results[0]
        resultsReads = results[1]
        resultsOnTarget = results[2]
        resultsCoverageLevels = results[3]
        resultsFirstColumn = results[4]
        resultsLastColumn = results[5]
        resultsUncountedStats = results[6]
        resultsInsertSizes = results[7]
        
        # Update uncounted stats
        self.uncounted["unmapped"] += resultsUncountedStats["unmapped"]
        self.uncounted["duplicate"] += resultsUncountedStats["duplicate"]
        self.uncounted["mapquality"] += resultsUncountedStats["mapquality"]
        
        # Update insert size stats, I stop counting after 10000000
        if len(self.insertSize) < 10000000:
            self.insertSize.extend(resultsInsertSizes)
        
        # We have to account for overlap of reads before adjusting total counts
        readOverlap = set(self.lastChunkColumn).intersection(set(resultsFirstColumn))
        if len(readOverlap) > 0:
            resultsReads -= len(readOverlap)
            for readId in readOverlap:
                for regionName in resultsFirstColumn[readId]:
                    # This takes care of situation where an overlapping read was counted as on-target for the same region set for both results
                    if regionName in self.lastChunkColumn[readId]:
                        resultsOnTarget[regionName] -= 1
        self.lastChunkColumn = resultsLastColumn
                        
        # With overlapping reads figured out, increment our totals
        self.totalReads += resultsReads
        for descriptor in resultsOnTarget.keys():
            self.onTarget[descriptor] += resultsOnTarget[descriptor]
        self.totalCoverage += resultsCoverageLevels[0]
        self.totalLength += resultsRegion.length
    
    def report(self, bamInput, genome = False):
        
        report = {}
        
        # Insert size calculation:    
        insertMean = sum(self.insertSize) / float(len(self.insertSize))
        insertSD = math.sqrt(sum(map(lambda x:(x-insertMean)**2, self.insertSize))/(len(self.insertSize)-1))
        
        allReads = self.totalReads
        
        for key,value in self.uncounted.items():
            allReads += value
                
        report["version"] = __version__
        report["inputBam"] = bamInput
        report["allReads"] = allReads
        report["readsCounted"] = self.totalReads
        report["readsNotCounted"] = self.uncounted
        report["insertMean"] = insertMean
        report["insertSD"] = insertSD
        report["onTarget"] = self.onTarget
        
        if genome:
            report["genome"] = { "avgCoverage" : float(self.totalCoverage) / self.totalLength }
        
        return report
    
    
    def __init__(self, regionSets):
        self.onTarget = {}
        if len(regionSets) > 0:
            for descriptor in regionSets:
                self.onTarget[descriptor] = 0.0
        
        self.totalReads = 0
        self.totalCoverage = 0
        self.totalLength = 0
        self.lastChunkColumn = {}
        self.uncounted = {"unmapped" : 0, "duplicate": 0, "mapquality": 0}
        self.insertSize = []          
            
    
    
class ProcessingRegionGenerator(object):
    
    def _sort(self, ):
        sortedRegions = {}
        
        def getStart(region):
            return region.start
        
        for chrom in self.regionByChromosome:
            sortedRegions[chrom] = sorted(self.regionByChromosome[chrom], key = getStart)
        
        self.regionByChromosome = sortedRegions
        self.sorted = True   
    
    def addRegion(self, region):
        region.index = self.regionCount
        self.regionCount += 1
        if region.chrom in self.regionByChromosome:
            self.regionByChromosome[region.chrom].append(region)
        else:
            self.regionByChromosome[region.chrom] = [region]
        
        self.sorted = False
    
    def returnProcessingRegion(self, ):
        '''Returns a list of tuples in the form [(region, [subRegion1, subRegion2...]), ...] where region and SubregionX are coveragekit.utils.region.Region objects.
        This groups regions specified by user into the processing chuncks definied by windowSize and genome size.
        
        :param header: pysam.AlignmentFile.header
        :type header: dict
        :param selectRegions: list of tuples returned by coveragekit.utils.bed.parseRegionBeds
        :type selectRegions: list
        :param windowSize: Size of bam chunk to be considered by a bam reader
        :type windowSize: int
        
        :rtype: dict
        
        '''
        if not self.sorted:
            self._sort()
        
        regionCount = 0
        processingRegions = []
        for sq in self.header['SQ']:
            chromLength = sq["LN"]
            chromName = sq["SN"]
            
            # The following deals with the stupid "chr" problem
            if chromName.startswith("chr"):
                editChromName = chromName[3:]
            else:
                editChromName = chromName
            
            chromStart = 0
            lastStop = 0
            if editChromName in self.regionByChromosome.keys():
                selectList = self.regionByChromosome[editChromName]
            else:
                selectList = []
    
            while lastStop < chromLength:
                subSelectRegions = []
                chromStart = lastStop
                lastStop += self.windowSize
                if lastStop > chromLength:
                    lastStop = chromLength
    
                curProcessingRegion = regionkit.Region(chromName, chromStart, lastStop, regionCount, "_processing", regionCount)
                regionCount += 1
                
                # Generate subregions by selections input by user
                selectCount = 0
                if (len(selectList) > 0) and (selectList[0].start < lastStop):
                    while (selectCount < len(selectList)):
                        if selectList[selectCount].stop < chromStart:
                            selectList.pop(selectCount)
                            continue
                        elif selectList[selectCount].start > lastStop:
                            break
                        elif selectList[selectCount].stop <= lastStop:
                            curSelect = selectList.pop(selectCount)
                        elif selectList[selectCount].stop > lastStop:
                            curSelect = selectList[selectCount]
                            selectCount += 1
    
                        if curSelect.start < chromStart:
                            selectStart = chromStart
                        else:
                            selectStart = curSelect.start
    
                        if curSelect.stop > lastStop:
                            selectStop = lastStop
                        else:
                            selectStop = curSelect.stop
    
                        subSelectRegions.append(regionkit.Region(editChromName, selectStart, selectStop, curSelect.name, curSelect.regionSet, curSelect.index))
    
                yield (curProcessingRegion, subSelectRegions)
    
    def __init__(self, bamFile, windowSize):
        bam = pysam.AlignmentFile(bamFile, 'rb')
        self.header = bam.header
        bam.close()
        self.regionByChromosome = {}
        self.windowSize = windowSize
        self.sorted = False
        self.regionCount = 0
    

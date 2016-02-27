
import json, logging

from coveragekit.version import __version__

class RegionSet(object):
    
    def add(self, region, levelReport):
        self.calcDone = False
        newRegion = True
        
        # Take care of pseudoautosomal
        if region.name in self.regionDict:
            if self.regionDict[region.name]["chrom"] == region.chrom:
                newRegion = False
            else:
                self.logger.warning("Potential ambiguity in gene name for {}. Chromosome {} versus {}.".format(region.name,region.chrom,self.regionDict[region.name]["chrom"]))
            
        if not newRegion:
            if region.regionSet != self.setName:
                raise Exception("Discordance between region set ({}) and {} of RegionSet object".format(region,self.setName))
            
            if region.start < self.regionDict[region.name]["start"]:
                self.regionDict[region.name]["start"] = region.start
            if region.stop > self.regionDict[region.name]["stop"]:
                self.regionDict[region.name]["stop"] = region.stop
            self.regionDict[region.name]["length"] += region.length
            self.regionDict[region.name]["coverage"] += levelReport[0]
            for curLevel in levelReport[1]:
                self.regionDict[region.name]["bg"][curLevel[2]].append(curLevel[:2])
            self.regionDict[region.name]["subregions"].append((region.start, region.stop, levelReport[0]))
        else:
            if region.regionSet != self.setName:
                raise Exception("Discordance between region set ({}) and {} of RegionSet object".format(region,self.setName))
            self.regionDict[region.name] = {}
            self.regionDict[region.name]["chrom"] = region.chrom
            self.regionDict[region.name]["start"] = region.start
            self.regionDict[region.name]["stop"] = region.stop
            self.regionDict[region.name]["length"] = region.length
            self.regionDict[region.name]["name"] = region.name
            self.regionDict[region.name]["index"] = self.numRegions
            self.regionDict[region.name]["coverage"] = levelReport[0]
            self.regionDict[region.name]["subregions"] = [(region.start,region.stop, levelReport[0])]
            self.regionDict[region.name]["bg"] = {}
            for i in self.levels:
                self.regionDict[region.name]["bg"][i] = []
            for curLevel in levelReport[1]:
                self.regionDict[region.name]["bg"][curLevel[2]].append(curLevel[:2])
            self.numRegions += 1
            
        self.coverage += levelReport[0]
        self.length += region.length
    
    def calc(self, ):
        for curRegionID,curRegion in self.regionDict.items():
            levelAggregate = 0
            curRegion["averageCoverage"] = curRegion["coverage"] / float(curRegion["length"])
            curRegion["levelCoverage"] = {}
            for i in reversed(self.levels):
                for j in curRegion["bg"][i]:
                    levelAggregate += j[1] - j[0]
                curRegion["levelCoverage"][i] = levelAggregate / float(curRegion["length"])
                self.levelCoverage[i] += levelAggregate
            #if levelAggregate < curRegion["length"]:
            #    print curRegion
            #    raw_input("Press ENTER to proceed")
        self.calcDone = True
    
    def report(self, ):
        if not self.calcDone:
            self.calc()
        report = {"name":self.setName,
                  "numRegions": self.numRegions,
                  "length": self.length,
                  "avgCoverage": self.coverage / float(self.length),
                  "coverageLevels": {}}
        for i in self.levels:
            report["coverageLevels"][i] = self.levelCoverage[i] / float(self.length)
        return report
    
    def _retrieve(self, regionID):
        if not self.calcDone:
            self.calc()
        try:
            r = self.regionDict[regionID]
        except KeyError:
            return None
        record = [regionID, r["chrom"], r["start"], r["stop"], json.dumps(r["subregions"]),r["length"],r["averageCoverage"], json.dumps(r["bg"])]
        for i in self.levels:
            record.append(r["levelCoverage"][i])
        record = tuple(record)
        return record
    
    def retrieve(self, regionID = None, regionList=None):
        if not self.calcDone:
            self.calc()
        if regionList is not None:
            retrieveList = regionList
        elif regionID is not None:
            retrieveList = [regionID]
        else:
            retrieveList = sorted(self.regionDict.keys())
        
        for r in retrieveList:
            yield self._retrieve(r)    
    
    def __init__(self, setName, levels):
        self.setName = setName
        self.length = 0
        self.coverage = 0
        self.numRegions = 0
        self.levels = tuple(sorted(levels))
        
        if self.levels[0] != 0:
            self.levels = (0,) + self.levels

        self.levelCoverage = {}
        for i in self.levels:
            self.levelCoverage[i] = 0
        
        self.regionDict = {}
        self.calcDone = False
        
        self.logger = logging.getLogger("coveragekit utils.region.RegionSet {}".format(self.setName))
        self.logger.setLevel(logging.INFO)
        

class Region(object):
    
    def __init__(self, chrom, start, stop, name, regionSet, index):
        self.chrom = str(chrom)
        self.start = int(start)
        self.stop = int(stop)
        self.name = str(name)
        self.regionSet = str(regionSet)
        self.index = int(index)
        
        self.length = stop - start
    
    def __repr__(self, ):
        return "{},{},{},{},{},{}".format(self.chrom,self.start,self.stop,self.name,self.regionSet,self.index)
    
#!/usr/bin/env python

import sqlite3, os, json, sys, optparse, logging

import coveragekit.utils.db
from coveragekit.utils.bed import stitchRegions
from coveragekit.version import __version__

def _prettifyResult(result, levelsMinKeys, levelsMaxKeys, dbLevels, reportRegions = True):
    prettifiedResult = {}
    prettifiedResult["id"] = result[0]
    prettifiedResult["position"] = "{}:{}-{}".format(result[1],result[2],result[3])
    prettifiedResult["coverage"] = result[6]
    
    prettifiedResult["percentGreaterOrEqual"] = {}
    for h in zip(dbLevels,result[8:]):
        prettifiedResult["percentGreaterOrEqual"][h[0]] = h[1]
    
    
    if reportRegions:
        prettifiedResult["coverageRegions"] = {"lessThan": {}, "greaterOrEqual": {}}
        resultLevels = json.loads(result[7])
        if len(levelsMaxKeys) > 0:
            for level in levelsMaxKeys:
                prettifiedResult["coverageRegions"]["lessThan"][level] = []
                toStitch = [] 
                for i in reversed(sorted(resultLevels.keys())):
                    if int(i) < level:
                        toStitch.extend(resultLevels[i])
                stitched = stitchRegions(toStitch)
                for j in stitched:
                    prettifiedResult["coverageRegions"]["lessThan"][level].append("{}:{}-{}".format(result[1],j[0],j[1]))
                    
        if len(levelsMinKeys) > 0:
            for level in levelsMinKeys:
                prettifiedResult["coverageRegions"]["greaterOrEqual"][level] = []
                toStitch = [] 
                for i in sorted(resultLevels.keys()):
                    if int(i) >= level:
                        toStitch.extend(resultLevels[i])
                stitched = stitchRegions(toStitch)
                for j in stitched:
                    prettifiedResult["coverageRegions"]["greaterOrEqual"][level].append("{}:{}-{}".format(result[1],j[0],j[1]))
    
    return prettifiedResult
    
def db(dbInput, genes = None, levelsMin = None, levelsMax = None, coverageMin = None, coverageMax = None, reportRegions = True):
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger("coveragekit db")
    logger.setLevel(logging.INFO)
    
    coverageDB = coveragekit.utils.db.CoverageDB(dbInput)
    
    if levelsMax:
        levelsMaxKeys = levelsMax.keys()
    else:
        levelsMaxKeys = []
    
    if levelsMin:
        levelsMinKeys = levelsMin.keys()
    else:
        levelsMinKeys = []
    
    dbLevelSet = set(coverageDB.levels)

    if (not set(levelsMaxKeys).issubset(dbLevelSet)) or (not set(levelsMinKeys).issubset(dbLevelSet)):
        print "The specified db only has the following levels available: {}".format(coverageDB.levels)
        sys.exit(1)
    
    returnedGenes = []
   
    results = {"meta" : {}, "queryResults" : [] }
    results["meta"]["version"] = __version__
    results["meta"]["dbSource"] = dbInput
    results["meta"]["coverageSource"] = coverageDB.coveragesource
    results["meta"]["regionSource"] = coverageDB.regionsource
    results["meta"]["dbLevels"] = coverageDB.levels
    for result in coverageDB.query(genes, coverageMin, coverageMax, levelsMin, levelsMax):
        results["queryResults"].append(_prettifyResult(result, levelsMinKeys, levelsMaxKeys, coverageDB.levels, reportRegions = reportRegions))
        returnedGenes.append(results["queryResults"][-1]["id"])
    
    results["meta"]["queryString"] = coverageDB.mostRecentQuery
    results["meta"]["queryResultNum"] = len(returnedGenes)
    
    if genes:
        notFound = set(genes)
        notFound.difference_update(set(returnedGenes))
        if notFound:
            logger.warning("The following regions were not found in the coverage database: [{}]".format(",".join(notFound)))
            
    results["queryResults"].sort(key=lambda k: k['id'])
    return results
    
def report(results, reportRegions = True, jsonOut = None, tsvOut = None):

    if jsonOut:
        with open(jsonOut, "w") as jsonFile:
            jsonFile.write(json.dumps(results, indent = 4, sort_keys = True))
            
    if tsvOut:
        with open(tsvOut, "w") as tsvFile:
            atOrAbove = []
            for l in results["meta"]["dbLevels"]:
                atOrAbove.append("PercentAtOrAbove{}X".format(l))
            
            tsvFile.write("RegionID\tPosition\tAverageCoverage\t{}".format("\t".join(atOrAbove)))
            if reportRegions:
                tsvFile.write("\tRegionsLessThan\tRegionsGreaterThanOrEqual\n")
            else:
                tsvFile.write("\n")
            for r in results["queryResults"]:
                atOrAbove = []
                for l in results["meta"]["dbLevels"]:
                    atOrAbove.append(str(r["percentGreaterOrEqual"][l]))
                tsvFile.write("{}\t{}\t{}\t{}".format(r["id"],r["position"],r["coverage"],"\t".join(atOrAbove)))
                
                if reportRegions:
                    tsvFile.write("\t{}\t{}\n".format(json.dumps(r["coverageRegions"]["lessThan"]), json.dumps(r["coverageRegions"]["greaterOrEqual"])))
                else:
                    tsvFile.write("\n")
            
    print "\n\ncoveragekit db results:"
    print "--------------"
    print "DB coverage source:\t{}".format(results["meta"]["coverageSource"])
    print "DB region source:\t{}".format(results["meta"]["regionSource"])
    print "DB query string:\t{}".format(results["meta"]["queryString"])
    print "Records retrieved:\t{}".format(results["meta"]["queryResultNum"])
    if jsonOut:
        print "JSON output:\t{}".format(jsonOut)
    if tsvOut:
        print "tsv output:\t{}".format(tsvOut)
    print "\n\n"

def run(inputArgs):
    usage = "%prog --db sample.db [ options ]"
    parser = optparse.OptionParser(usage=usage, prog = "coveragekit db")
    parser.add_option("-d","--db", dest="db", help="Input database.", default="")
    parser.add_option("--geneList", type="string", dest="geneList", help="Comma-separated gene list.", default="")
    parser.add_option("--geneListFile", type="string", dest="geneListFile", help="File with newline-separated gene list.", default="")
    parser.add_option("--levelsMin", type="string", dest="levelsMin", help="Comma-separated list of minimum percents at X coverage with colon delimited coverage level prepended ( eg '5:99,10:95,20:90' ).", default="")
    parser.add_option("--levelsMax", type="string", dest="levelsMax", help="Comma-separated list of maximum percents at X coverage with colon delimited coverage level prepended ( eg '5:99,10:95,20:90' ).", default="")
    parser.add_option("--coverageMin", type="float", dest="coverageMin", help="Minimum average coverage.", default=None)
    parser.add_option("--coverageMax", type="float", dest="coverageMax", help="Maximum average coverage.", default=None)
    parser.add_option("--reportRegions", action="store_true", dest="reportRegions", help="Report regions with coverage of intersest as JSON stings.", default=False)
    parser.add_option("--json", type="string", dest="json", help="Output JSON file.", default=None)
    parser.add_option("--tsv","--txt", type="string", dest="tsv", help="Output tsv file.", default=None)
    (options, args) = parser.parse_args(inputArgs)

    # Bam file is required as well as one output
    if len(options.db) == 0: parser.error("Missing db, use -d or --db.")
    if (options.json is None) and (options.tsv is None): parser.error("Must specify an output with --json or --tsv")
    
    if (len(options.geneList) > 0) and (len(options.geneListFile) > 0):
        parser.error("Cannot specify both --geneList and --geneListFile.")
    elif len(options.geneList) > 0:
        genes = options.geneList.split(",")
    elif len(options.geneListFile) > 0:
        genes = []
        with open(options.geneListFile) as geneFile:
            for line in geneFile:
                genes.append(line.rstrip())
    else:
        genes = None
    
    if len(options.levelsMin) > 0:
        levelsMin = {}
        for level in options.levelsMin.split(','):
            key,value = level.split(':',1)
            levelsMin[int(key)] = float(value)
    else:
        levelsMin = None
    
    if len(options.levelsMax) > 0:
        levelsMax = {}
        for level in options.levelsMax.split(','):
            #print level
            key,value = level.split(':',1)
            levelsMax[int(key)] = float(value)
    else:
        levelsMax = None
    
    results = db(options.db, genes = genes, levelsMin = levelsMin, levelsMax = levelsMax, coverageMin = options.coverageMin, coverageMax = options.coverageMax, reportRegions = options.reportRegions)
    report(results, reportRegions = options.reportRegions, jsonOut = options.json, tsvOut = options.tsv)

if __name__ == '__main__':
    run(sys.argv[1:])
from coveragekit.version import __version__
from coveragekit.utils.region import Region

def bedToRegions(descriptor, bedFile):
    with open(bedFile) as regionFile:
        regionCount = 0
        for line in regionFile:
            # Get rid of godforsaken "chr" that some insist on adding to chromosome name
            if line.startswith("chr"):
                line = line[3:]

            lineSplit = line.rstrip().split("\t")
            try:
                regionName = lineSplit[3]
            except IndexError:
                regionName = regionCount

            newRegion = Region(lineSplit[0], int(lineSplit[1]), int(lineSplit[2]), regionName, descriptor, regionCount)
            regionCount += 1
            yield newRegion

def stitchRegions(unstitched):
    unstitched.sort(key=lambda x: int(x[0]))
    stitched = []
    if len(unstitched) > 0:
        lastInterval = unstitched.pop(0)
        lastStart = lastInterval[0]
        lastStop = lastInterval[1]
            
        for i in unstitched:
            if i[0] == lastStop:
                lastStop = i[1]
            else:
                stitched.append([lastStart, lastStop])
                lastStart = i[0]
                lastStop = i[1]
        stitched.append([lastStart, lastStop])
    return stitched
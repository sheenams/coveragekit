
from coveragekit.version import __version__

class CoverageLevel(object):
    
    def _add(self, pos, coverage):
        self.coverage += coverage
        newLevel = self.curLevel
        while (coverage >= self.curLevelMax):
            newLevel += 1
            self.curLevelMin = self.levels[newLevel]
            self.curLevelMax = self.levels[newLevel + 1]
        while (coverage < self.curLevelMin):
            newLevel -= 1
            self.curLevelMin = self.levels[newLevel]
            self.curLevelMax = self.levels[newLevel + 1]
        if newLevel != self.curLevel:
            if self.curPos >= self.curLevelStart:
                coverageRegion = (self.curLevelStart, pos)
                self.coverageRegions[self.levels[self.curLevel]].append(coverageRegion)
            self.curLevelStart = pos
            self.curLevel = newLevel
        self.curPos = pos
               
    def add(self, pos, coverage):
        if pos > (self.curPos + 1):
            for i in range(self.curPos + 1, pos):
                self._add(i, 0)
        elif pos < (self.curPos + 1):
            #print pos, self.curPos, self.start, self.stop
            raise Exception("CoverageLevel.add can only go left to right along chromosome.")
        self._add(pos, coverage)
            
    def report(self, ):        
        # Close out currently open region
        coverageRegion = (self.curLevelStart, self.stop)
        self.coverageRegions[self.levels[self.curLevel]].append(coverageRegion)

        bg = []
        for l in reversed(self.levels[:-1]): # We don't want to include the inf coverage levels
            for c in self.coverageRegions[l]:
                '''if c[0] < self.start:
                    start = self.start
                else:
                    start = c[0]'''
                bg.append((c[0],c[1],l))
        bg.sort(key=lambda x: x[0])
        return (self.coverage, bg)    
    
    def __init__(self, start, stop, levels):
        self.start = start
        self.stop = stop
        self.levels = levels
        self.coverage = 0
        self.coverageRegions = {}
        
        if self.levels[0] != 0:
            self.levels = (0,) + self.levels
        if self.levels[-1] != float("Inf"):
            self.levels = self.levels + (float("Inf"),)
            
        for i in self.levels[:-1]:
            self.coverageRegions[i] = []
        
            
        self.curLevelMin = self.levels[0]
        self.curLevelMax = self.levels[1]
        self.curLevelStart = start
        self.curLevel = 0
        self.curPos = start - 1
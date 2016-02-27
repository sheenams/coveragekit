
import sqlite3, os, json, logging, datetime

from coveragekit.version import __version__

class CoverageDB(object):
    
    def _create(self, dbfile, regionsource, coveragesource, levels, mapq, dup, totalCoverage):
        self.conn = sqlite3.connect(dbfile)
        self.regionsource = regionsource
        self.coveragesource = coveragesource
        self.levels = levels
        self.c = self.conn.cursor()
        
        # Create basic table for regions (genes)
        self.c.execute("CREATE TABLE regions(id text, chrom text, start integer, stop integer, subregions text, length integer, coverage real, levels text)")
        
        # Add additional database columns depending on the coverage levels assayed
        for i in self.levels:
            self.c.execute("ALTER TABLE regions ADD COLUMN 'percent{}X' real".format(i))
            #self.c.execute("ALTER TABLE regions ADD COLUMN 'level{}regions' text".format(i))
        
        # Index on region (gene) id
        self.c.execute("CREATE UNIQUE INDEX ididx ON regions(id)")
        self.c.execute("CREATE INDEX coverageidx ON regions(coverage)")
        
        # Create metadata table and put in the appropriate data
        self.c.execute("CREATE TABLE metadata(regionsource text, coveragesource text, levels text, mapqualityCutoff int, duplicatesAllowed int, totalCoverage int)")
        levelString = ",".join(str(x) for x in levels)
        self.c.execute("INSERT INTO metadata VALUES (?,?,?,?,?,?)",(regionsource, coveragesource, levelString, mapq, dup, totalCoverage))
        
        # Create coveragekit metadata table and put in the appropriate data
        self.c.execute("CREATE TABLE coveragekit(version text, dateCreated text)")
        self.c.execute("INSERT INTO coveragekit VALUES (?,?)",(__version__, datetime.datetime.now().isoformat()))
        
        self.conn.commit()
    
    def _load(self, dbfile):
        self.conn = sqlite3.connect(dbfile)
        self.c = self.conn.cursor()
        self.c.execute("SELECT * FROM metadata LIMIT 1")
        metadata = self.c.fetchone()
        self.regionsource = metadata[0]
        self.coveragesource = metadata[1]
        self.levels = tuple([int(x) for x in metadata[2].split(",")])
        
    def insert(self, region, levels):
        pass
    
    def insertRegionSet(self, regionSet):
        for setRecord in regionSet.retrieve():
            columns = ",".join("?"*(8+(len(self.levels))))
            self.c.execute("INSERT INTO regions VALUES ({})".format(columns),setRecord)
        self.conn.commit()
    
    def query(self, geneID = None, coverageLowCutoff = None, coverageHighCutoff = None, levelsLowCutoff = None, levelsHighCutoff = None):
        queryString = []
        if geneID:
            queryString.append('id IN ("{}")'.format('","'.join(geneID)))
        
        if coverageLowCutoff:
            queryString.append("coverage >= {}".format(coverageLowCutoff))
            
        if coverageHighCutoff:
            queryString.append("coverage < {}".format(coverageHighCutoff))
            
        if levelsLowCutoff:
            levelQuery = []
            for curLevel,cutoff in levelsLowCutoff.items():
                if cutoff != '.':
                    levelQuery.append("percent{}X >= {}".format(curLevel,cutoff/100.0))
            queryString.append(" AND ".join(levelQuery))
        
        if levelsHighCutoff:
            levelQuery = []
            for curLevel,cutoff in levelsHighCutoff.items():
                if cutoff != '.':
                    levelQuery.append("percent{}X < {}".format(curLevel,cutoff/100.0))
            queryString.append(" AND ".join(levelQuery))
        self.logger.debug('SELECT * FROM regions WHERE {}'.format(" AND ".join(queryString)))
        self.mostRecentQuery = 'SELECT * FROM regions WHERE {}'.format(" AND ".join(queryString))
        return self.c.execute('SELECT * FROM regions WHERE {}'.format(" AND ".join(queryString)))            
    
    def reset(self, regionsource, coveragesource, levels, mapq, dup, totalCoverage):
        self.conn.close()
        self.regionsource = None
        self.coveragesource = None
        self.levels = None
        os.remove(self.dbfile)
        self._create(self.dbfile, regionsource, coveragesource, levels, mapq, dup, totalCoverage)    
    
    def __init__(self, db, regionsource = None, coveragesource = None, levels = None, mapq = 1, dups = False, totalCoverage = 0, overwrite = False):
        self.dbfile = db
        if os.path.isfile(db):
            self._load(db)
            if overwrite:
                if regionsource and coveragesource and levels:
                    self.reset(regionsource, coveragesource, levels, mapq, dups, totalCoverage)
                else:
                    raise Exception("Cannot create coverage database without regionsource, bamsource and levels specified.")
        else:
            if regionsource and coveragesource and levels:
                self._create(db, regionsource, coveragesource, levels, mapq, dups, totalCoverage)
            else:
                raise Exception("Cannot create coverage database without regionsource, bamsource and levels specified.")
        self.logger = logging.getLogger("coveragekit utils.db.CoverageDB")
        self.logger.setLevel(logging.INFO)
        self.mostRecentQuery = ""
    
    def __del__(self, ):
        self.conn.close()
    
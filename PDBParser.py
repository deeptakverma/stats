import sys, math

class PDBParse:
    def __init__(self, File):
        self.Coord = []
        ResNumList =[]
        AAConverter = {"ALA": "A","GLU": "E","GLN": "Q","ASP": "D","ASN": "N","LEU": "L","GLY": "G","LYS": "K","SER": "S","VAL": "V","ARG": "R","THR": "T","PRO": "P","ILE": "I","MET": "M","PHE": "F","TYR": "Y","CYS": "C","CYX": "C","CYD": "C","TRP": "W","HID": "H","HIE": "H","HIS": "H","HIP":"H"}
        self.seq = ""
        for i in open(File).readlines():
            if i[:4] == "ATOM":
                l = LineParse(i)
                ResCheck = False
                if l.ResName == "GLY":
                    if l.AtomName == "CA":
                        self.Coord.append(l)
                        self.seq += AAConverter[l.ResName]
                        ResNumList.append(l.ResNum)
                        ResCheck = True
                else:
                    if l.AtomName == "CB":
                        self.Coord.append(l)
                        self.seq += AAConverter[l.ResName]
                        ResNumList.append(l.ResNum)
                        ResCheck = True
                if len(ResNumList) > 1 and ResCheck:
                    if ResNumList[-1]-ResNumList[-2] != 1:
                        sys.stderr.write("This protein is not continuous. Check at %d\n"%ResNumList[-1])
                        sys.exit()
        self.Length = len(self.seq)
        self.StartingResidue = ResNumList[0]

        cutoff = 8
        self.ContactMap = []
        for i in xrange(self.Length):
            if self.seq[i] == "P":
                CM_line = [0]*self.Length
            else:
                CM_line = []
                for j in xrange(self.Length):
                    if self.seq[j] == "P":
                        CM_line.append(0)
                    else:
                        a = self.Coord[i]; b= self.Coord[j]
                        distance = math.sqrt((a.X-b.X)**2+(a.Y-b.Y)**2+(a.Z-b.Z)**2)
                        if distance < cutoff: CM_line.append(1)
                        else: CM_line.append(0)
            self.ContactMap.append(CM_line)

        self.Pairs = []
        for i in xrange(self.Length-1):
            for j in xrange(i,self.Length):
                if i==j: pass
                else:
                    if self.ContactMap[i][j]: self.Pairs.append((i,j))

class LineParse:
    def __init__(self, line):
        self.AtomName = line[12:16].strip()
        self.ResName = line[17:20].strip()
        self.ResNum = int(line[22:26])
        self.X = float(line[30:38])
        self.Y = float(line[38:46])
        self.Z = float(line[46:54])

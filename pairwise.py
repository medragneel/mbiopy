from proteins_masses import  proteins_codes
import re

class My_seq:
    def __init__(self,seq,biotype="DNA"):
        self.seq = seq.upper()
        self.biotype = biotype
    def __len__(self):
        return len(self.seq)
    def __getitem__(self,n):
        return self.seq[n]
    def __getslice__(self,i,j):
        return self.seq[i:j]
    def __str__(self):
        return self.seq
    def get_biotype(self):
        return self.biotype
    def __info__(self):
        return f"sequence: {self.seq} \n biotype: {self.biotype}"
    def alphabet (self):
      if (self.biotype=="DNA"): return "ACGT"
      elif (self.biotype=="RNA"): return "ACGU"
      elif (self.biotype=="PROTEIN"): return "ACDEFGHIKLMNPQRSTVWY"
      else: return ""

    def validate (self):
        alp = self.alphabet()
        res = True
        i = 0
        while i < len(self.seq) and res:
            if self.seq[i] not in alp: res = False
            else: i += 1
        return res 
    def transcribe_dna(self):
        if self.biotype == "DNA":
            return My_seq(self.seq.replace("T","U"), "RNA")
        else:
            return None
    def reverse_comp (self):
        if (self.biotype != "DNA"): return None
        comp = ""
        for c in self.seq:
            if (c == 'A'): comp = "T" + comp 
            elif (c == "T"): comp = "A" + comp 
            elif (c == "G"): comp = "C" + comp
            elif (c== "C"): comp = "G" + comp
        return My_seq(comp, "DNA")

    def translate_dna(self):
        peptide=[]
        triplets=re.findall('...?',self.seq)
        for n in triplets:
            for key,_ in proteins_codes.items():
                if n in proteins_codes[key]:
                    peptide.append(key)
        aa = ''.join(peptide)   

        return My_seq(aa,"PROTEIN")







class Align:
    def __init__(self,ls,al_type="P"):
        self.ls = ls
        self.al_type = al_type 
        self.alphabet = ""
        self.sm = {}
    def __len__(self):
        return len(self.ls[0])
    def __getitem__(self,n):
        if type(n) is tuple and len(n) == 2:
            i,j = n
            return self.ls[i][j]
        elif type(n) is int: return self.ls[n]
        return None
    def __str__(self):
        res=""
        for sq in self.ls:
            res += "\n" + sq
        return res
    def num_seqs(self):
        return len(self.ls)
    def col(self,idx):
        return [self.ls[k][idx] for k in range(len(self.ls))]
    def consensus(self):
        cons=""
        for i in range(len(self)):
            cont={}
            for k in range(len(self.ls)):
                c = self.ls[k][i]
                if c in cont:
                    cont[c] = cont[c] + 1
                else:
                    cont[c] = 1
            maximum=0
            cmax=None
            for key in cont.keys():
                if key != "−" and cont[key] > maximum:
                    maximum= cont[key]
                    cmax= key
            cons += str(cmax)
        return cons



class Submat:
    def __init__(self):
        self.alphabet = ""
        self.sm = {}
    def __getitem__(self,ij):
        i,j=ij
        return self.score_pair(i,j)

     
    def score_pair(self,c1,c2):
        if c1 not in self.alphabet or c2 not in self.alphabet:
            return None
        return self.sm[c1+c2]
        
     
    def create_sub_mat(self,m,mm,alphabet):
        self.alphabet = alphabet
        for c1 in alphabet:
            for c2 in alphabet:
                if (c1==c2):
                    self.sm[c1+c2] = m
                else:
                   self.sm[c1+c2] = mm
        return None

    def read_sub_mat_from_file(self,filename,sep):
        with open(filename ,"r") as blosum62:
            line=blosum62.readline().split(sep)
            self.alphabet = ""
            for i in range(0,len(line)):
                self.alphabet += line[i][0]
            
            for i in range(0,len(line)):    
                ln = blosum62.readline().split(sep)
                for j in range(0,len(ln)):
                    k = self.alphabet[i] + self.alphabet[j]
                    self.sm[k] = int(ln[j])
            return None



class Pairwise:
    def __init__(self,sm,g):
        self.sm = sm 
        self.g = g
        self.S = {}
        self.T = {}
        self.s1 = ""
        self.s2 = ""
     
    def score_pos(self,c1,c2):
        if c1  == "−" or c2 == "−":
            return self.g
        else:
            return self.sm[c1+c2]
    def score_alin (self, al):
        res = 0;
        for i in range(len(al)):
            res += self.score_pos (al[0][i], al[1][i])
        return res
    def nw(self,s1,s2):
        if (s1.biotype != s2.biotype): return None
        self.S = [[0]]
        self.T = [[0]]
        self.s1 = s1
        self.s2 = s2
        for j in range(1,len(s2) + 1):
            self.S[0].append(self.g * j)
            self.T[0].append(3)
        for i in range(1,len(s1) + 1):
            self.S.append([self.g * i])
            self.T.append([2])
        for i in range(0,len(s1)):
            for j in range(len(s2)):
                rs1= self.S[i][j] + self.score_pos(s1[i],s2[j])
                rs2 = self.S[i][j+1] + self.g
                rs3 = self.S[i+1][j] + self.g
                self.S[i+1].append(max(rs1,rs2,rs3))
                self.T[i+1].append(max3t(rs1,rs2,rs3))
        return self.S[len(s1)][len(s2)]
    def final_align(self):
        res=["",""]
        i = len(self.s1)
        j= len(self.s2)
        while i > 0 or j > 0:
            if self.T[i][j] == 1:
                res[0] =  self.s1[i-1] + res[0]
                res[1] = self.s2[j-1] + res[1]
            elif self.T[i][j] == 3:
                res[0] = "-" + res[0]
                res[1] = self.s2[j-1] + res[1] 
                j -= 1
            elif self.T[i][j] ==2:
                res[0] = self.s1[i-1] + res[0]
                res[1] = "-" + res[1]
                i -= 1
            else:break
        return Align(res,self.s1.biotype)
    def sw(self,s1,s2):
        if (s1.biotype != s2.biotype): return None
        self.S = [[0]]
        self.T = [[0]]
        self.s1 = s1
        self.s2 = s2
        maxscore = 0 
        # column insirtion
        for _ in range(1,len(s2) + 1 ):
            self.S[0].append(0)
            self.T[0].append(0)
        #row shit
        for _ in range(1,len(s1) + 1):
            self.S.append([ 0 ])
            self.T.append([ 0 ])
        for i in range(0,len(s1)):
            for j in range(0,len(s2)): 
                rs1 = self.S[i][j] + self.score_pos(s1[i],s2[j])
                rs2 = self.S[i][j+1] + self.g
                rs3 = self.S[i+1][j] + self.g
                b= max(rs1,rs2,rs3)
                if(b <= 0):
                    self.S[i+1].append(0)
                    self.T[i+1].append(0)
                else:
                    self.S[i+1].append(b)
                    self.T[i+1].append(max3t(rs1,rs2,rs3))
                    if b > maxscore:
                        maxscore = b
        return maxscore



    def recover_align_local(self):
        res = ["",""]
        maxscore = 0
        maxrow = 0
        maxcol = 0
        for i in range(1,len(self.S)):
            for j in range(1,len(self.S[i])):
                if self.S[i][j] > maxscore:
                    maxscore = self.S[i][j]
                    maxrow = i
                    maxcol = j
        i = maxrow
        j=maxcol
        while i > 0 or j > 0:
            if self.T[i][j]==1:
                res[0] = self.s1[i-1] + res[0]
                res[1] = self.s2[j-1] +res[1]
                i -= 1
                j-= 1
            elif self.T[i][j] == 3 :
                res[0] = "-" + res[0]
                res[1] = self.s2[j-1] +res[1]
                j-=1
            elif self.T[i][j] == 2:
                res[0] = self.s1[i-1] +res[0]
                res[1] = "-" + res[1]
                i -=1
            else: break
        return Align(res,self.s1.biotype)



def max3t (v1, v2, v3):
    if v1 > v2:
        if v1 > v3: return 1
        else: return 3
    else:
        if v2 > v3: return 2
        else: return 3


def print_mat(mat):
    for i in range (0,len(mat)):
        print(mat[i])

def test_seq():
    s1 = My_seq("ATGTGATAAGAATAGAATGCTGAATAAATAGAATGACAT")
    s2 = My_seq("MKVVLSVQERSVVSLL", "PROTEIN")
    print(s1.validate(), s2.validate())
    print(s1.transcribe_dna().translate_dna())

def test_align():
    al= Align(ls=["ATGA−A","AA−AT−"],al_type="D")
    print(len(al))
    print (al.__len__())
    print(al.col(idx=1))
    print(al[1,1])
    print(al.consensus())

def test_pairwise():
    seq1 = My_seq("ATGATATGATGATT")
    seq2 = My_seq("GATGAATAGATGTGT")
    sm = Submat()
    sm.create_sub_mat(3, -1, "ACGT")
    alin = Pairwise(sm, -3)
    print("------- smith and waterman------------")
    print(alin.sw(seq1, seq2))
    print_mat(alin.S)
    print(alin.recover_align_local())
    
    print("------- needleman and wunch -------------")
    print(alin.nw(seq1,seq2))
    print_mat(alin.S)
    print(alin.final_align())

test_seq()
test_align()
test_pairwise()

import re
from proteins_masses import pms,proteins_codes


class Btp:
    def Gc(self,seq):
        counter_gc=0
        for n in seq.upper():
            if n == 'C' or n =='G':
                counter_gc+=1 

        return (counter_gc / len(seq)) * 100
        
    def transcribe_dna(self,seq):
        return re.sub('T','U',seq.upper())
    def dna_reverse(self,seq):
        reverse=[]
        for n in seq.upper():
            if n == 'A':
                n = 'T'
                reverse.append(n)
            elif n == 'T':
                n = 'A'
                reverse.append(n)
            elif n == 'C':
                n = 'G'
                reverse.append(n)

            elif n == 'G':
                n = 'C'
                reverse.append(n)
        #print(seq)
        #print('|' * len(seq))


        return ''.join(reverse)
    def find_mut(self,seq1,seq2):
        list_mut=[]
        count = 0
        len_sq= 0
        if len(seq1) < len(seq2):
            len_sq=len(seq1)
        elif len(seq1) > len(seq2):
            len_sq=len(seq2)
        elif len(seq1) == len(seq2):
            len_sq= len(seq1)
        for i in range(len_sq):
            if seq1[i].lower() != seq2[i].lower():
                count+=1
                print(i)
                print(seq1[i],seq2[i])
                t_mut=(i,seq1[i],seq2[i])
                list_mut.append(t_mut)
        return list_mut
    def length(self,seq):
        return len(seq)
    def calculate_pms(self,seq):
        pmass=0
        for key,value in pms.items():
            for p in seq.upper():
                if key == p:
                    pmass += value

        return pmass
    def translate_dna(self,seq):
        peptide=[]
        triplets=re.findall('...?',self.transcribe_dna(seq.replace(' ','')))
        print(triplets)
        for n in triplets:
            for key,_ in proteins_codes.items():
                if n in proteins_codes[key]:
                    print(n)
                    peptide.append(key)
            

        print(peptide)
        return''.join(peptide)



               

        

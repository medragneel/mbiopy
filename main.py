from flask import Flask,render_template,request
from btp import Btp



app= Flask(__name__,template_folder='templates')

bio=Btp()
gc_per = 0
rna_seq=''
cmp= ''
cmp_dna_seq=''
seq1=''
seq2=''
mut_count=0
sq= ''
len_sq= 0
ps=''
pm=0
psq=''
pep=''
mut_len= 0

@app.route('/')
def home():
    return render_template('index.html')


@app.route('/tools')
def tools():
    return render_template('tools.html')

@app.route('/tools/gc',methods=['GET','POST'])
def gc():
    global gc_per
    if request.method == 'POST':
        gc_seq= request.form['gc_seq']
        gc_per= bio.Gc(gc_seq) 
    return render_template('gc.html',per= gc_per)

@app.route('/tools/rna',methods=['GET','POST'])
def rna():
    global rna_seq
    if request.method == 'POST':
        dna_seq= request.form['dna']
        rna_seq= bio.transcribe_dna(dna_seq) 
    return render_template('rna.html',rna=rna_seq)

@app.route('/tools/dna_complementary',methods=['GET','POST'])
def dna_cmp():
    global cmp,cmp_dna_seq
    if request.method == 'POST':
        cmp_dna_seq= request.form['dna_cmp']
        cmp= bio.dna_reverse(cmp_dna_seq) 
    return render_template('dna_comp.html',f=cmp_dna_seq,cmp=cmp)

@app.route('/tools/mutations_count',methods=['GET','POST'])
def mutations_count():
    global seq1,seq2,mut_count,mut_len
    if request.method == 'POST':
        seq1=request.form['seq1']
        seq2=request.form['seq2']
        print(seq1)
        print(seq2)
        mut_count= bio.find_mut(seq1,seq2)
        mut_len= len(mut_count)
    return render_template('mutation.html',seq1=seq1,seq2=seq2,mc=mut_count,l=mut_len)

@app.route('/tools/seq_length',methods=['GET','POST'])
def seq_length():
    global sq,len_sq
    if request.method == 'POST':
        sq=request.form['len']
        len_sq = bio.length(sq)
    return render_template('length.html',sq=sq,length=len_sq)

@app.route('/tools/pmass',methods=['GET','POST'])
def pmass():
    global ps,pm
    if request.method == 'POST':
        ps=request.form['pmass']
        pm = bio.calculate_pms(ps)
    return render_template('pmass.html',ps=ps,pm=pm)

@app.route('/tools/protein_seq',methods=['GET','POST'])
def pseq():
    global psq,pep
    if request.method == 'POST':
        psq=request.form['protein_seq']
        pep= bio.translate_dna(psq)
    return render_template('prot_seq.html',psq=psq,pep=pep)
@app.route('/about')
def about():
    return render_template('about.html')



if __name__ == "__main__":
    app.run(debug=True)





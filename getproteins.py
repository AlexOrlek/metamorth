import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from mymod import runsubprocess

#arg[1] is outputpath; arg[2] is genbank file

#read through each plasmid sequence record, extract proteins, write protein fasta for each plasmid
outputpath=sys.argv[1]
genbankfile=sys.argv[2]

runsubprocess(['mkdir -p %s/fastas'%outputpath],shell=True)

f3=open('%s/fastas/sequencestats.tsv'%outputpath,'w')
f3.write('%s\t%s\t%s\t%s\t%s\t%s\n'%('accession','num_pseudogenes','num_genes','num_genes_excludingskipped','num_genes_missingid','num_genes_missingtranslation'))
f4=open('%s/fastafilepaths.tsv'%outputpath,'w')
f5=open('%s/blastdbfilepaths.tsv'%outputpath,'w')

with open(genbankfile) as f:
    for indx, seq_record in enumerate(SeqIO.parse(f, "genbank")):
        header=str(seq_record.id)
        print header
        f2=open('%s/fastas/%s.fasta'%(filepath,header),'w')
        pseudocounter=0
        genecounter=0
        skippedgenemissingidcounter=0
        skippedgenemissingtranslationcounter=0
        #if indx>1:
        #    break
        #print header
        if seq_record.features:
            for feature in seq_record.features:
                if feature.type=="CDS":
                    if 'pseudo' in feature.qualifiers or 'pseudogene' in feature.qualifiers: #'product' not in
                        print('gene is pseudogene')
                        pseudocounter=pseudocounter+1
                        continue
                    genecounter=genecounter+1
                    #if counter>10:
                    #    break
                    try:
                        protid=feature.qualifiers["protein_id"]
                    except:
                        try:
                            print('protein id missing for accession: %s locus tag: %s; skipping'%(header, feature.qualifiers["locus_tag"]))
                        except:
                            print('protein id missing for accession: %s'%header)
                        skippedgenemissingidcounter=skippedgenemissingidcounter+1
                        continue
                    assert len(protid)==1, "there are multiple protein ids associated with feature"
                    protid=protid[0]
                    try:
                        genename=feature.qualifiers['gene'][0]
                    except:
                        genename='-'
                    try:
                        protname=feature.qualifiers['product'][0]
                    except:
                        protname='-'

                    strand=feature.strand
                    start=int(feature.location.start)
                    end=int(feature.location.end)                    
                    #protnuclseq=feature.location.extract(seq_record).seq
                    #try:
                    #    protseq1=protnuclseq.translate(table=11,cds=True)
                    #except:
                    #    print('translation error; skipping')
                    #    continue

                    try:
                        protseq2=feature.qualifiers['translation']
                    except:
                        print('no translation for %s|%s; skipping'%(header, protid))
                        skippedgenemissingtranslationcounter=skippedgenemissingtranslationcounter+1
                        continue
                    assert len(protseq2)==1, "there are multiple translations associated with feature"
                    protseq2=Seq(protseq2[0])
                    #assert protseq1==protseq2,"there is a disagreement between translated nucleotide sequence and provided translated protein sequence"
                    print '%s|%s protein:%s gene:%s'%(header,protid, protname, genename)
                    #print protseq2
                    
                    #make sequence record and write to file
                    seqid='%s|%s|%s|%s-%s'%(header,protid,strand,start,end)
                    seqdescr='protein:%s gene:%s'%(protname,genename)
                    seq_record2=SeqRecord(protseq2)
                    seq_record2.id=seqid
                    seq_record2.description=seqdescr
                    SeqIO.write(seq_record2,f2,"fasta")
        f2.close()
        skippedgenecounter=skippedgenemissingidcounter+skippedgenemissingtranslationcounter
        f3.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(header,pseudocounter,genecounter,genecounter-skippedgenecounter,skippedgenemissingidcounter,skippedgenemissingtranslationcounter))
        fastafilepath='%s/fastas/%s.fasta'%(filepath,header)
        f4.write('%s\t%s\n'%(header,fastafilepath))
        blastdbfilepath='%s/fastas/%s_db'%(filepath,header)
        f5.write('%s\t%s\n'%(header,blastdbfilepath))
f3.close()
f4.close()
f5.close()

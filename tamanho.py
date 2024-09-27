from Bio import SeqIO
import pandas as pd


# 1. Ler a tabela de classificação original
classification_df = pd.read_csv('/home/hppp123/IC/tamanho/classificationtable_OG_TESTE_tab.tsv', sep='\s+', header=None)
classification_df.columns = ['Group', 'Orthogroup', 'ID']


#print(classification_df.head())
#print(f"Número de colunas : {classification_df.shape[1]}")


# 2. dicionários para armazenar o comprimento dos transcritos e CDS
transcript_lengths = {}
cds_lengths = {}

# 3. Calcular o comprimento dos transcritos
with open('/home/hppp123/IC/tamanho/ORTHO_TESTE.fasta') as transcripts_fasta:
    for record in SeqIO.parse(transcripts_fasta, "fasta"):
        transcript_lengths[record.id] = len(record.seq)

# 4. Calcular o comprimento das CDS
with open('/home/hppp123/IC/tamanho/CDS_TESTE.fasta') as cds_fasta:
    for record in SeqIO.parse(cds_fasta, "fasta"):
        cds_lengths[record.id] = len(record.seq)

# 5. DataFrame com Group, Orthogroup, ID, Transcript_Length e CDS_Length
new_df = classification_df[['Group', 'Orthogroup', 'ID']].copy()  # Copiar as 3 primeiras colunas
new_df['Transcript_Length'] = new_df['ID'].map(transcript_lengths)  # Adicionar a coluna de comprimentos dos transcritos
new_df['CDS_Length'] = new_df['ID'].map(cds_lengths)  # Adicionar a coluna de comprimentos dos CDS

# 6. tabela com as colunas criadas
new_df.to_csv('Classification_lengths.tsv', sep='\t', index=False)

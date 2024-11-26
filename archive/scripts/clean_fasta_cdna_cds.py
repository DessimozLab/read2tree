
from Bio import SeqIO
from Bio.Seq import Seq
import sys
from os import listdir
import os




def read_fasta_files(input_folder, format_input="fna"):

    files = listdir(input_folder)
    records_all = []
    file_names = [] 
    for file in files:
        sp_name = file.split(".")[:-1]
        if file.split(".")[-1] == format_input:
            file_names.append(file)
            records = list(SeqIO.parse(input_folder + file, "fasta"))
            records_all.append(records)        
        else:
            print("we are not reading the file "+str(input_folder+file)+" since extension is not faa.")
    if records_all:
        print("there are ", len(file_names), format_input, " files, and the first file has ", len(records_all[0]), "sequences in it.") 
    else:
        print("there is no  " +format_input, " files in ",input_folder) 
    return file_names, records_all


def create_five_letter(file_names, output_five_letter_tsv = "clean_five_letter_species.tsv"):
    
    fiveLetter_species_dic = {}
    countr = 0
    for file_name in file_names:
        fiveLetter_species = "s" + str(countr).zfill(4) 
        fiveLetter_species_dic[file_name] = fiveLetter_species
        countr += 1
    file_out = open(output_five_letter_tsv, "w")
    for species_name, fiveLetter in fiveLetter_species_dic.items():
        file_out.write(species_name + "\t" + fiveLetter + "\n")
    file_out.close()
    print("the five letter codes for each faa files are written in "+output_five_letter_tsv)

    return fiveLetter_species_dic



def clean_translate(records ,species_fivelet):
    
    records_nuc = []
    records_aa = []
    for record in records:
        sequence = record.seq
        remainder = len(sequence) % 3
        if remainder != 0:
            sequence +=Seq('N' * (3 - remainder)) 
            record.seq= sequence
    
        id_old = str(record.id).replace("_","").replace(".","")
        id_new=  species_fivelet + id_old
        
        nuc_seq= SeqIO.SeqRecord(sequence, id=id_new, description="cleaned for r2t", name = id_new)
        
        protein_seq = sequence.translate()
        protein_seq = SeqIO.SeqRecord(protein_seq, id=id_new, description="cleaned for r2t", name = id_new)
    
        
        records_nuc.append(nuc_seq)
        records_aa.append(protein_seq)
    
    print("the clean aa and nuc for "+species_fivelet+" is ready")
    
    return records_nuc, records_aa





if __name__ == '__main__':

    input_folder_fna = sys.argv[1] + "/"  # "myfolder/input_fna/" #
        
    file_names, records_all = read_fasta_files(input_folder_fna, "fna")
    fiveLetter_species_dic = create_five_letter(file_names)
    
    
    folder_aa= "clean_aa"
    
    
    if not os.path.exists(folder_aa):
        os.makedirs(folder_aa)
    else:
        print("ERROR the folder exists "+folder_aa +" better to remove it ")
    
    records_nuc_all_clean=[]
    for idx in range(len(file_names)):
        file_name = file_names[idx]
        records = records_all[idx]    
        species_fivelet = fiveLetter_species_dic[file_name]
    
        records_nuc, records_aa = clean_translate(records ,species_fivelet)
            
        SeqIO.write(records_aa, folder_aa+"/"+species_fivelet+".fa", "fasta")
        
        records_nuc_all_clean += records_nuc # one big list 
    
    
    SeqIO.write(records_nuc_all_clean, "dna_ref.fa", "fasta")
    
    print("we wrote "+str(len(file_names))+" faa fiels in the folder "+folder_aa+" and the nucluetide sequences all together in dna_ref.fa" )

    print("Now you can use the folder with OMA standalone" )

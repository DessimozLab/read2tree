
from Bio import SeqIO
import sys
from os import  listdir

def read_fasta_files(input_folder_faa,format_input):
    files = listdir(input_folder_faa)
    fa_all = []
    species_name_all = [] # which are the file name of fasta files

    for file in files:
        sp_name = file.split(".")[:-1]
        if file.split(".")[-1] == format_input:  
            species_name_all.append(".".join(sp_name))
            records_prot = list(SeqIO.parse(input_folder_faa+file, "fasta"))
            fa_all.append(records_prot)
    print("there are ",len(species_name_all),format_input, " files, and the first file has ",len(fa_all[0]),"sequences in it.")     # , sum([len(i) for i in og_all]

    return (species_name_all, fa_all)

def read_fiveLetter_species_file(input_five_letter_csv):
    fiveLetter_species_dic ={}
    file1 = open(input_five_letter_csv,"r")
    for line in file1:
        species_name, fiveLetter_species = line.strip().split("\t")
        fiveLetter_species_dic[species_name] = fiveLetter_species
    file1.close()
    return fiveLetter_species_dic


def write_fiveLetter_species_file(species_name_all, output_five_letter_tsv):
    
    fiveLetter_species_dic = {}
    try:
        for species_name in species_name_all: # let's try to extract a code  which unique from file name
            fiveLetter_species= species_name.split(".")[0].split("_")[1][-5:]  # GCA_000849305.1_ViralProj14697_translated_cds.faa  # JN032115.1_cds_from_genomic.fna
        fiveLetter_species_dic[species_name] = fiveLetter_species
    except:
        fiveLetter_species_dic = {}

    if len(set(fiveLetter_species_dic.values())) != len(set(species_name_all)):
        #"we assume the last five letter of NCBI is unique, please provide five letter code for species name as input as  Read2tree works with five letter code specicies name." 
        fiveLetter_species_dic = {}
        countr=0
        for species_name in species_name_all:
            fiveLetter_species= "s"+str(countr).zfill(4) #species_name.split(".")[0].split("_")[1][-5:]
            fiveLetter_species_dic[species_name] = fiveLetter_species
            countr+=1
    
    
    file1 = open(output_five_letter_tsv,"w")
    for species_name, fiveLetter in fiveLetter_species_dic.items():
        file1.write(species_name+"\t"+fiveLetter+"\n")    
    file1.close()
    
    return fiveLetter_species_dic




def edit_record_write_fna(species_name_all_fna, fna_all, output_file_fna):
    # Add the five letter species code to each record in fasta file
    all_prot_fna = []
    for species_idx_fna, species_record in enumerate(fna_all):
        species_name_fna = species_name_all_fna[species_idx_fna]
        species_name_faa = species_name_fna[:-17]+"_translated_cds"
        fiveLetter_species = fiveLetter_species_dic[species_name_faa]
        #print(species_name_fna, fiveLetter_species)

        for prot in species_record: 
                                     
            prot_id_old = prot.id.split(" ")[0][3:]
            # lcl|KY249672.1_prot_APW78783.1_1 [protein=NS1] [protein_id=APW78783.1] [location=99..518] [gbkey=CDS] 
            # >lcl|AF092942.1_cds_AAC96311.1_11
            # >lcl|AF092942.1_prot_AAC96311.1_11 
            prot_id_old_split= prot_id_old.split("_")
            try:
                prot_id_old_split.remove("prot")
            except:
                try:
                    prot_id_old_split.remove("cds")
                except:
                    print("Error: prot/cds is not inside the record id of ", prot.id.split(" ")[0])
                    print("we expect such format >lcl|AF092942.1_cds_AAC96311.1_11. Contact the developers.")
                    exit            
            
            prot_id_edit = ".".join(prot_id_old_split)
            
            prot_id_new = fiveLetter_species+ prot_id_edit
            prot.id = prot_id_new
            prot.name = prot_id_new
            prot.description = prot_id_new


            all_prot_fna.append(prot)

    SeqIO.write(all_prot_fna, output_file_fna, "fasta")
    
    return all_prot_fna


def edit_record_write_faa(species_name_all_faa, faa_all, fiveLetter_species_dic, output_folder_faa,all_prot_fna_id_set):

    # Add the five letter species code to each record in fasta file
    for species_idx, species_record in enumerate(faa_all):
        species_name = species_name_all_faa[species_idx]
        fiveLetter_species = fiveLetter_species_dic[species_name]
        for prot in species_record:    
            prot_id_old = prot.id.split(" ")[0][3:]
            # lcl|KY249672.1_prot_APW78783.1_1 [protein=NS1] [protein_id=APW78783.1] [location=99..518] [gbkey=CDS] 
            # >lcl|AF092942.1_cds_AAC96311.1_11
            # >lcl|AF092942.1_prot_AAC96311.1_11 
            prot_id_old_split= prot_id_old.split("_")
            try:
                prot_id_old_split.remove("prot")
            except:
                try:
                    prot_id_old_split.remove("cds")
                except:
                    print("Error: prot/cds is not inside the record id of ", prot.id.split(" ")[0])
                    print("we expect such format >lcl|AF092942.1_cds_AAC96311.1_11. Contact the developers. ")
                    exit
            prot_id_edit = ".".join(prot_id_old_split)

            prot_id_new = fiveLetter_species+ prot_id_edit
            
            assert prot_id_new in all_prot_fna_id_set,prot_id_old+"is not in fna file (exact match after removing _cds_ or _prot_)"
            prot.id = prot_id_new
            prot.name = prot_id_new
            prot.description = prot_id_new


        SeqIO.write(species_record, output_folder_faa+fiveLetter_species+".fa", "fasta")
        
    return faa_all
    



if __name__ == '__main__':

    input_folder_faa = sys.argv[1]+"/"  # "data/" 
    output_folder_faa = sys.argv[2]+"/" # "DB/" 

    output_file_fna = sys.argv[3]   # "all_cdna.fa" 
    
    if len(sys.argv)>4:
        input_five_letter_tsv = sys.argv[4]
    else:
        input_five_letter_tsv = ""    
    
    
    '''
    $ cat five_letter_species.tsv 
    GCA_003266525.1_ASM326652v1_translated_cds	66525
    GCA_000857565.1_ViralProj15251_translated_cds	57565
    GCA_000849305.1_ViralProj14697_translated_cds	49305
    
    '''
    
    output_five_letter_tsv =  input_folder_faa+"five_letter_species.tsv" # argv[2]


    input_folder_fna = input_folder_faa

    (species_name_all_faa, faa_all)  = read_fasta_files(input_folder_faa,"faa")
    (species_name_all_fna, fna_all)  = read_fasta_files(input_folder_fna,"fna")

    assert len(species_name_all_faa) ==len(species_name_all_fna), "the number of faa and fna files should be the same in the folder."

    assert len(faa_all) ==len(fna_all), "the number of faa and fna records should be the same."

    if input_five_letter_tsv:
        fiveLetter_species_dic = read_fiveLetter_species_file(input_five_letter_tsv)

    else:
        fiveLetter_species_dic=  write_fiveLetter_species_file(species_name_all_faa, output_five_letter_tsv)



    all_prot_fna = edit_record_write_fna(species_name_all_fna, fna_all, output_file_fna)
    print("Edited cdna records are written to the file",output_file_fna)

    all_prot_fna_recordid = [i.id for i in all_prot_fna]
    all_prot_fna_id_set = set(all_prot_fna_recordid)
    assert len(all_prot_fna_recordid) == len(all_prot_fna_id_set), "all record id in fna files should be unique. we consider this format when we checl"  +all_prot_fna_recordid[0]


    faa_all = edit_record_write_faa(species_name_all_faa, faa_all, fiveLetter_species_dic, output_folder_faa, all_prot_fna_id_set)

    print("Edited protien records are written to the folder",output_folder_faa)




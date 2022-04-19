from tables import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyoma.browser import db
import familyanalyzer as fa

# parameters
MIN_SPECIES = 20
DUP_RATIO = 0
DIR = '/Users/daviddylus/Research/read2tree/reference_datasets/Dataset1/Output/'

# read in files
hog_XML = DIR+'HierarchicalGroups.orthoxml'
og_XML = DIR+'OrthologousGroups.orthoxml'
h5file = open_file("/Volumes/Untitled/OmaServer.h5", mode="r")

genomeTab = h5file.root.Genome
dbObj = db.Database(h5file)
omaIdObj = db.OmaIdMapper(dbObj)

if DUP_RATIO != 0:
  hog_op = fa.OrthoXMLParser(hog_XML)
  gene_family_xml_nodes_hog = hog_op.getToplevelGroups()
  # select all the families with more than X species and duplication ratio smaller than Y
  hog_families_X = {}
  for i, family in enumerate(gene_family_xml_nodes_hog):
    family_id = family.get('id')
    genes_per_hog = [val for sublist in hog_op.getGenesPerSpeciesInFam(family).values() for val in sublist]
    species_per_hog = hog_op.getGenesPerSpeciesInFam(family).keys()
    duplication_ratio = float(len(genes_per_hog)) / float(len(species_per_hog))
    if len(species_per_hog) >= MIN_SPECIES and duplication_ratio <= DUP_RATIO:
      hog_families_X[family_id] = genes_per_hog

  print(len(hog_families_X))


og_op = fa.OrthoXMLParser(og_XML)
gene_family_xml_nodes_og = og_op.getToplevelGroups()
og_families_X = {}
for i, family in enumerate(gene_family_xml_nodes_og):
    family_id = family.get('id')
    genes_per_og = [val for sublist in og_op.getGenesPerSpeciesInFam(family).values() for val in sublist]
    species_per_og = og_op.getGenesPerSpeciesInFam(family).keys()
    if len(species_per_og) >= MIN_SPECIES:
      og_families_X[family_id] = genes_per_og

print(len(og_families_X))

if DUP_RATIO != 0:
  family_map = {}
  entries_map_omaids = {}
  cpt = 0
  for og in og_families_X:
    cpt += 1
    if cpt % 10 == 0:
      print("{} on {}".format(cpt, len(og_families_X)))
    a = og_families_X[og]
    for hog in hog_families_X:
      b = hog_families_X[hog]
      if len(set(a).intersection(b)) == 30:
        oma_ids_full = [og_op.mapGeneToXRef(val, 'protId') for val in og_families_X[og]]
        oma_ids = [og_op.mapGeneToXRef(val, 'protId').split(' | ')[0] for val in og_families_X[og]]
        entries = [omaIdObj.omaid_to_entry_nr(val) for val in oma_ids]
        for oma_id in oma_ids_full:
          entries_map_omaids[omaIdObj.omaid_to_entry_nr(oma_id.split(' | ')[0])] = oma_id
        family_map[og] = entries
        break
  print(len(entries_map_omaids))
else:
  family_map = {}
  entries_map_omaids = {}
  cpt = 0
  for og in og_families_X:
    cpt += 1
    if cpt % 1000 == 0:
      print(og_op.mapGeneToXRef(og_families_X[og][0], 'protId').split(' | ')[0])
      print("{} on {}".format(cpt, len(og_families_X)))
    oma_ids_full = [og_op.mapGeneToXRef(val, 'protId') for val in og_families_X[og] if og_op.mapGeneToXRef(val, 'protId')]
    oma_ids = [val.split(' | ')[0] for val in oma_ids_full]
    entries = [omaIdObj.omaid_to_entry_nr(val) for val in oma_ids if omaIdObj.omaid_to_entry_nr(val)]
    print(entries)
    for oma_id in oma_ids_full:
      entries_map_omaids[omaIdObj.omaid_to_entry_nr(oma_id.split(' | ')[0])] = oma_id
    family_map[og] = entries
  print(len(entries_map_omaids))


family_map_invert = {}
for key in family_map:
  for val in family_map[key]:
    family_map_invert[val]=key

print(len(family_map_invert))

records = []
for key in family_map_invert:
  new_id = entries_map_omaids[key] + '| OG' + family_map_invert[key]
  record = SeqRecord(Seq(dbObj.get_cdna(key)), id=new_id, description="")
  records.append(record)

with open("dataset2.fasta", "w") as output_handle:
  SeqIO.write(records, output_handle, "fasta")


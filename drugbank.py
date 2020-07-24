# To run:
# PYTHONIOENCODING=UTF-8 python3 drugbank.py

import os
import re
import requests
import wget
import xml.etree.ElementTree as ET
from datetime import date
from atomwrappers import *

xml_file = "raw_data/drugbank/full database.xml"
tag_prefix = "{http://www.drugbank.ca}"
output_file = "dataset/drugbank_{}.scm".format(str(date.today()))

xml_root = ET.parse(xml_file).getroot()

if os.path.exists(os.path.join(os.getcwd(), output_file)):
  os.remove(output_file)
out_fp = open(output_file, "a", encoding = "utf8")

def find_tag(obj, tag):
  return obj.find(tag_prefix + tag)

def findall_tag(obj, tag):
  return obj.findall(tag_prefix + tag)

def get_child_tag_text(obj, tag):
  return find_tag(obj, tag).text

def find_mol_type(mol):
  if "CHEBI:" in mol.upper():
    mol_type = ChebiNode(mol)
  elif "PubChem:" in mol or "PubChemSID" in mol:
    mol_type = PubchemNode(mol)
  else:
    mol_type = CMoleculeNode(mol)
  return mol_type

def get_pubchem_cid(sid):
  print("--- Getting PubChem CID for SID:{}\n".format(sid))
  try:
    response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/" + sid + "/cids/txt", timeout=20)
  except:
    print("=== Connection error")
    return None

  if response.status_code != 200:
    print("=== Failed to find a PubChem CID for SID:{}\n".format(sid))
    return None
  else:
    return response.text.strip()

# Get ChEBI IDs for reference later
chebi_obo = "raw_data/chebi.obo"
chebi_url = "ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo"

if os.path.exists(chebi_obo):
  print("Removing file: {}".format(chebi_obo))
  os.remove(chebi_obo)

chebi_file = wget.download(chebi_url, "raw_data")
print("\nFile downloaded: {}".format(chebi_file))

chebi_fp = open(chebi_file, "r", errors="ignore")
chebi_dict = {}
chebi_name = []
chebi_id = None

for line in chebi_fp:
  line = line.replace("\n", "")
  if line == "[Term]":
    if len(chebi_name) > 0 and chebi_id != None:
      for name in chebi_name:
        chebi_dict[name.lower()] = "ChEBI:" + chebi_id
      # print("ChEBI ID: {}\nName: {}\n".format(chebi_id, chebi_name))
    chebi_name = []
    chebi_id = None
  elif line.startswith("id: "):
    chebi_id = line.replace("id: CHEBI:", "")
  elif line.startswith("name: "):
    chebi_name.append(line.replace("name: ", ""))
  elif line.startswith("synonym: ") and "EXACT" in line:
    name = re.match(".+\"(.+)\".+", line).group(1)
    if name not in chebi_name:
      chebi_name.append(name)

chebi_fp.close()

# Then go through the whole file once, to get the external IDs
id_dict = {}
for drug in xml_root:
  drugbank_id = get_child_tag_text(drug, "drugbank-id")
  name = get_child_tag_text(drug, "name").lower()

  chebi = None
  pubchem_cid = None
  pubchem_sid = None

  for external_id in findall_tag(find_tag(drug, "external-identifiers"), "external-identifier"):
    resource = get_child_tag_text(external_id, "resource")
    identifier = get_child_tag_text(external_id, "identifier")
    if resource == "ChEBI":
      chebi = "ChEBI:" + identifier
    elif resource == "PubChem Compound":
      pubchem_cid = "PubChem:" + identifier
    elif resource == "PubChem Substance":
      # Prefix will be added later
      pubchem_sid = identifier

  # Try to get the ChEBI ID from the official database if it's not found in DrugBank
  if chebi == None:
    chebi = chebi_dict.get(name)

  # Try to get the PubChem CID from the official database if it's not found in DrugBank
  if pubchem_cid == None and pubchem_sid != None:
    pubchem_cid = get_pubchem_cid(pubchem_sid)

  if chebi != None:
    id_dict[drugbank_id] = chebi
  elif pubchem_cid != None:
    id_dict[drugbank_id] = pubchem_cid
  elif pubchem_sid != None:
    id_dict[drugbank_id] = "PubChemSID:" + pubchem_sid
  else:
    # If no desired external IDs is found, use the DrugBank ID
    id_dict[drugbank_id] = "DrugBank:" + drugbank_id

# Finally do the conversion for each of the drugs
drug_groups = []
for drug_tag in xml_root:
  drugbank_id = get_child_tag_text(drug_tag, "drugbank-id")
  standard_id = id_dict.get(drugbank_id)
  name = get_child_tag_text(drug_tag, "name").lower()
  description = get_child_tag_text(drug_tag, "description")

  standard_id = find_mol_type(standard_id)
  evalink = CEvaluationLink(CPredicateNode("has_name"), CListLink(standard_id, CConceptNode(name)))
  out_fp.write(evalink.recursive_print() + "\n")

  if description != None:
    description = description.replace("\"", "\\\"").strip()
    evalink = CEvaluationLink(CPredicateNode("has_description"), CListLink(standard_id, CConceptNode(description)))
    out_fp.write(evalink.recursive_print() + "\n")

  for group_tag in findall_tag(find_tag(drug_tag, "groups"), "group"):
    drug_group = group_tag.text + " drug"
    inhlink = CInheritanceLink(standard_id, CConceptNode(drug_group))
    out_fp.write(inhlink.recursive_print() + "\n")
    if drug_group not in drug_groups:
      inhlink = CInheritanceLink(CConceptNode(drug_group), CConceptNode("drug"))
      out_fp.write(inhlink.recursive_print() + "\n")
      drug_groups.append(drug_group)

  general_references_tag = find_tag(drug_tag, "general-references")
  articles_tag = find_tag(general_references_tag, "articles")
  for article_tag in findall_tag(articles_tag, "article"):
    pubmed_id = get_child_tag_text(article_tag, "pubmed-id")
    if pubmed_id != None:
      pubmed_id = "https://www.ncbi.nlm.nih.gov/pubmed/?term=" + pubmed_id
      evalink = CEvaluationLink(CPredicateNode("has_pubmedID"), CListLink(standard_id, CConceptNode(pubmed_id)))
      out_fp.write(evalink.recursive_print() + "\n")

  drug_interactions_tag = find_tag(drug_tag, "drug-interactions")
  for drug_interaction_tag in findall_tag(drug_interactions_tag, "drug-interaction"):
    other_drug_drugbank_id = get_child_tag_text(drug_interaction_tag, "drugbank-id")
    other_drug_standard_id = id_dict.get(other_drug_drugbank_id)
    # For some reason a few of them are not in the 'full database' file?
    if other_drug_standard_id == None:
      other_drug_standard_id = other_drug_drugbank_id
    
    other_drug_standard_id = find_mol_type(other_drug_standard_id)
    evalink = CEvaluationLink(CPredicateNode("interacts_with"), CListLink(standard_id, other_drug_standard_id))
    out_fp.write(evalink.recursive_print() + "\n")

  pathways_tag = find_tag(drug_tag, "pathways")
  for pathway_tag in findall_tag(pathways_tag, "pathway"):
    smpdb_id = get_child_tag_text(pathway_tag, "smpdb-id")
    for involved_drug_tag in findall_tag(find_tag(pathway_tag, "drugs"), "drug"):
      involved_drug_drugbank_id = get_child_tag_text(involved_drug_tag, "drugbank-id")
      involved_drug_standard_id = id_dict.get(involved_drug_drugbank_id)
      # For some reason a few of them are not in the 'full database' file?
      if involved_drug_standard_id == None:
        involved_drug_standard_id = involved_drug_drugbank_id
      
      involved_drug_standard_id = find_mol_type(involved_drug_standard_id)
      memberlink = CMemberLink(involved_drug_standard_id, SMPNode(smpdb_id))
      out_fp.write(memberlink.recursive_print() + "\n") 

    for uniprot_id_tag in findall_tag(find_tag(pathway_tag, "enzymes"), "uniprot-id"):
      uniprot_id = uniprot_id_tag.text
      evalink = CEvaluationLink(CPredicateNode("catalyzed_by"), CListLink(SMPNode(smpdb_id), ProteinNode(uniprot_id)))
      out_fp.write(evalink.recursive_print() + "\n")

  targets_tag = find_tag(drug_tag, "targets")
  for target_tag in findall_tag(targets_tag, "target"):
    be_id = get_child_tag_text(target_tag, "id")
    polupeptide_tag = find_tag(target_tag, "polypeptide")
    uniprot_id = polupeptide_tag.attrib["id"] if polupeptide_tag else None
    name = get_child_tag_text(target_tag, "name").strip().lower()
    action_tags = findall_tag(find_tag(target_tag, "actions"), "action")
    # Some drug has an unknown action yet not marked as "unknown", use "unknown" as well for them
    action = action_tags[0].text if action_tags else "unknown"
    target_id = "Uniprot:" + uniprot_id if uniprot_id else "DrugBank:" + be_id

    # TODO: Generate as directional (ListLink) for all of them for now
    target_id = find_mol_type(target_id)
    evalink = CEvaluationLink(CPredicateNode(action), CListLink(standard_id, target_id))
    out_fp.write(evalink.recursive_print() + "\n")
    evalink = CEvaluationLink(CPredicateNode("has_name"), CListLink(target_id, CConceptNode(name)))
    out_fp.write(evalink.recursive_print() + "\n")

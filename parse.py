import xml.etree.ElementTree as ET
import pandas as pd

def get_atc(drug):
    atc_codes = drug.findall(r'./{http://www.drugbank.ca}atc-codes//{http://www.drugbank.ca}atc-code')
    atc_codes = [elem.get("code") for elem in atc_codes if elem.get("code") is not None]
    return '|'.join(atc_codes)

def elements2string(target, elements):
    return '\n'.join([elem.text for elem in elements.findall(target) if elem.text is not None])

def small_molecule(result_path=r'N_ATC.xlsx', drugbank_path=r"./full database.xml", small_molecule=True, intersection=False):
    """Extract information from drugbank's xml file.

    Args:
        result_path (str, optional): File path for result excel. Defaults to r'N_ATC.xlsx'.
        drugbank_path (str, optional): Drugbank's full database .xml file. Defaults to r"./full database.xml".
        small_molecule (bool, optional): Only parse small molecule drugs. Defaults to True.
        intersection (bool, optional): whether to get intersection. Defaults to False.
    """
    
    print('Start parsing...')
    tree = ET.parse(drugbank_path)
    root = tree.getroot()
    print('All drugs', len(root))
    if small_molecule:
        drugs = root.findall(r'./{http://www.drugbank.ca}drug[@type="small molecule"]')
        print("Number of small molecule drugs: ", len(drugs))
    else: 
        drugs = root.findall(r'./{http://www.drugbank.ca}drug')
    
    approved_num=0
    code_num=0
    category_num=0
    intersection=0
    
    unii1=0
    unii2=0
    unii3=0
    unii4=0
    unii5=0
    
    target_df = pd.DataFrame(columns=['drugbank_id', 'name', 'atc_codes'])
    n_df = pd.DataFrame(columns=['drugbank_id', 'name', 'atc_codes'])
    for index, drug in enumerate(drugs):
        group = elements2string(r'.//{http://www.drugbank.ca}group', drug)
        unii = drug.findall(r'./{http://www.drugbank.ca}unii')
        if len(unii) == 1 and None not in [elem.text for elem in unii]:
            unii5+=1
        if 'approved' not in group:
            continue
        
        approved_num += 1
        unii = drug.findall(r'./{http://www.drugbank.ca}unii')
        if len(unii) == 1 and None not in [elem.text for elem in unii]:
            unii4+=1
        # N drug
        code = drug.findall(r'.//{http://www.drugbank.ca}level[@code="N"]')
        if code:
            code_num+=1
            unii = drug.findall(r'./{http://www.drugbank.ca}unii')
            n_df.loc[index, 'drugbank_id'] = elements2string(r'./{http://www.drugbank.ca}drugbank-id[@primary="true"]', drug)
            n_df.loc[index, 'name'] = elements2string(r'./{http://www.drugbank.ca}name', drug)
            n_df.loc[index, 'atc_codes'] = get_atc(drug)

            if len(unii) == 1 and None not in [elem.text for elem in unii]:
                unii1+=1
            
        # CENTRAL NERVOUS SYSTEM AGENTS
        category = elements2string(r'./{http://www.drugbank.ca}categories//{http://www.drugbank.ca}category', drug)
        if "CENTRAL NERVOUS SYSTEM AGENTS" in category.upper():
            category_num+=1
            unii = drug.findall(r'./{http://www.drugbank.ca}unii')
            if len(unii) == 1 and None not in [elem.text for elem in unii]:
                unii2+=1
            
        # Intersection
        if not intersection:
            continue 
        code = drug.findall(r'.//{http://www.drugbank.ca}level[@code="N"]')
        if code:
            category = elements2string(r'./{http://www.drugbank.ca}categories//{http://www.drugbank.ca}category', drug)
            if "CENTRAL NERVOUS SYSTEM AGENTS" in category.upper():
                intersection+=1
                target_df.loc[index, 'drugbank_id'] = elements2string(r'./{http://www.drugbank.ca}drugbank-id[@primary="true"]', drug)
                target_df.loc[index, 'name'] = elements2string(r'./{http://www.drugbank.ca}name', drug)
                target_df.loc[index, 'atc_codes'] = get_atc(drug)
                calculated_properties = drug.find(r'./{http://www.drugbank.ca}calculated-properties')
                for properties in calculated_properties:
                    target_df.loc[index, properties.find(r"{http://www.drugbank.ca}kind").text] = properties.find(r"{http://www.drugbank.ca}value").text

                unii = drug.findall(r'./{http://www.drugbank.ca}unii')
                if len(unii) == 1 and None not in [elem.text for elem in unii]:
                    unii3+=1

            
            
    print("Number of approved drugs: ", approved_num)
    print("Number of approved drugs with atc code N: ", code_num)
    print("Number of approved drugs with category: ", category_num)
    if intersection:
        print("Number of drugs in intersection: ", intersection)
    print("Uni:\n\tN ATC: {}, CENTRAL NERVOUS SYSTEM AGENTS:{}, Intersection:{}, Approved:{}, Small molecule drugs:{}".format(unii1, unii2, unii3, unii4, unii5))
    n_df.to_excel(r'N_ATC.xlsx', index=False, columns=["drugbank_id", "name", "atc_codes"])
            
    
if __name__ == "__main__":  
    small_molecule(
        result_path=r'N_ATC.xlsx', 
        drugbank_path=r"./full database.xml", 
        small_molecule=True, 
        intersection=False
        )
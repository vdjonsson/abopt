import structure 
import pandas as pd 


def test(inpt):
    mut_panda = pd.read_csv(inpt)
    mut_dict = {}
    for index, row in mut_panda.iterrows():
        if not row['wildtype']:
            pdb_name = row['pdb']
            mut = row['mut']
            if pdb_name in mut_dict:
                mut_dict[pdb_name].append(mut)
            else:
                mut_dict[pdb_name] = [mut]
    for pdb in mut_dict:
        temp = mut_panda.loc[mut_panda['pdb'] == pdb]
        temp_ab = temp.iloc[0, :]
        ab = temp_ab['antibody']
        indiv_list = data_path + 'individual_list_' + pdb + '.txt'
        write_to_mutations_file(mut_dict[pdb], indiv_list)
        pdb_url = 'https://files.rcsb.org/download/' + pdb + '.pdb'
        pdb_file = wget.download(pdb_url)
        run_repair_model(ab, pdb.lower())
        rep_pdb = pdb.lower() + '_Repair'
        run_build_model(ab, rep_pdb, indiv_list)
        rename_bm_out(rep_pdb, indiv_list)
        for mut in mut_dict[pdb]:
            mut_pdb = rep_pdb + '_' + mut
            print(mut_pdb)
            print(label_chains(pdb.lower()))
            rm_virus_with_biopandas(mut_pdb, label_chains(pdb.lower()))
            no_vir_pdb = mut_pdb + '_no_virus'
            run_repair_model(ab, no_vir_pdb)

# inp = '/home/teafs/Downloads/tmp.csv'
# test(inp)

import MLPopt_suppotfunc as MSF

pathtxt = 'cal/output/mol_to_ad/floder_name.txt'
pathads = 'cal/output/mol_to_ad/'
folder_dict = MSF.txt_to_dict(pathtxt)
namelist = list(folder_dict.keys())
smileslist =list(folder_dict.values())
for foldername in namelist:
    fp = pathads+foldername
    print(fp)
    files = MSF.find_files(fp,'.vasp')
    print(files)


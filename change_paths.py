import os
current_dir = os.getcwd()
with open("/home/user/cosmo/wCDM_training/lua_files/test_0.lua", "r") as f, \
     open("mod.lua", "w") as f_mod:
    lines = f.readlines()
    new_file_lines = lines
    for i, line in enumerate(lines):
        if "scratch" in line:
            new_file_lines[i] = line.replace("/scratch/decola/joao.reboucas2/COLA_projects/wCDM", current_dir)
    f_mod.writelines(new_file_lines)
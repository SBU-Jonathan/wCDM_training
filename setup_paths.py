"""
Description: Adjusts Santos Dumont paths in lua files and transferinfo files to your current directory.

Author: João Rebouças, 29-05-2023
"""
import os

current_dir = os.getcwd()

def setup_current_path_in_file(file_path):
    """
    Description: Changes the hard-coded string "/scratch/decola/joao.reboucas2/COLA_projects/wCDM" to the current directory.
    """
    with open(file_path, "r") as f:
        lines = f.readlines()
        new_file_lines = lines
        for i, line in enumerate(lines):
            if "scratch" in line: # Santos Dumont paths start with /scratch/...
                new_file_lines[i] = line.replace("/scratch/decola/joao.reboucas2/COLA_projects/wCDM", current_dir)
    
    with open(file_path, "w") as f:
         f.writelines(new_file_lines)

if __name__ == "__main__":
    print("Adjusting lua files and transferinfo files...")
    for i in range(500):
        setup_current_path_in_file(f"./transferinfo_files/transferinfo{i}.txt")
        setup_current_path_in_file(f"./lua_files/run_{i}.lua")
        setup_current_path_in_file(f"./lua_files/run_{i}ph_rev.lua")
        
    for i in range(3):
        setup_current_path_in_file(f"./transferinfo_files/transferinfotest_{i}.txt")
        setup_current_path_in_file(f"./lua_files/test_{i}.lua")
        setup_current_path_in_file(f"./lua_files/test_{i}ph_rev.lua")
    
    print("Finished! You may now run the simulations.")
    
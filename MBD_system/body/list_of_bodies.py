"""
Created on 18. feb. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import logging
import os
import re


def create_list(file_abs_path):
    """
    Lists file bodies.txt in folder - folder_path and returns the list of body data files.
    Args:
        filename - absolute path of folder
        folder_name - name of the folder
    Returns:
        
    """
    bodies = []
    if os.path.isfile(file_abs_path):
        with open(file_abs_path, "r") as file_:
        #    'r' - read
        #    'w' - write
        #    'a' - append
        #    'r+' - read and write
        
            for line in file_:
                if line.startswith('BODIES FILE-START'):
                    continue

                elif line.startswith('BODIES FILE-END'):
                    file_.close()
                    break

                elif line.startswith("#"):
                    pass

                else:
                    match_line_end_ = re.search(r'\n', line)
                    body_filename = line[0:match_line_end_.start()]
                    
                    bodies.append(os.path.splitext(body_filename)[0])
            
        # logging.getLogger("DyS_logger").info("File %s found! Bodies created successfully!", filename)
                    
    else:
        logging.getLogger("DyS_logger").info("File %s not found! No bodies created!", filename)
    
    return bodies

if __name__ == "__main__":
    print create_list("bodies.txt")
#     #    list all files in path = abs_path
#     list_of_files = os.listdir(MBD_folder_abs_path)
#     
#     #    file extension
#     extension = "stl"
#     
#     #    predefine list
#     list_of_names = []
#     for filename in list_of_files:
#         filename_ = filename.lower()
# 
#         
#         #    search where the extension in filename starts
#         match_extension = re.search(r".", filename_)
#         #    removes the extension
#         filename_without_extension = filename_[0:match_extension.start()-4]
#         
#         data_filename = "data_"+filename_without_extension+".txt"
# 
#         if filename_.endswith(extension) and os.path.isfile(os.path.join(MBD_folder_abs_path, data_filename)):
#             #    saves the filename without extension to list
#             list_of_names.append(filename_without_extension)
#     
#     return list_of_names


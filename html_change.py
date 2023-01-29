# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 06:19:55 2022

@author: User
"""
import pandas

####################################################
# Decide which files need to be updated
####################################################
# reading the CSV file
config = pandas.read_csv("change_config.csv")

names = config['dir'] + config['filename']
base_dir = config['base_dir']
titles = config['title']
update_files = config['update']

# names = ["examples/markets", "examples/operation", "examples/framework"]
# base_dir = ["../", "../", "../"]
# titles = ["Results: Power Market", "Results: System Operation", "Parameters"]

update_index = input("Update index.html? yes = 1; no = 0 | ")
print("\n")
# update_files = []
# for iter in range(len(names)):
#     temp = input("Update " + names[iter] + ".html? yes = 1; no = 0 | ")     
#     update_files.append(temp)


####################################################
# Code for html content alternation of index file
####################################################
if update_index == "1":
    filename = "index.html"
    title = "trEnD: Project Documentation"
    # Read in the file
    with open(filename, 'r') as file:
      filedata = file.read()
    
      
    # Find and elinimate the code for header
    find_string = "<div id=\"header\">\n\n\n\n"
    find_string += "<h1 class=\"title toc-ignore\">"
    find_string += title
    find_string += "</h1>\n\n</div>" 
    
    replace_string = ""
    filedata = filedata.replace(find_string, replace_string)

    # Write the file out 
    with open(filename, 'w') as file:
        file.write(filedata)
        

####################################################
# Code for html content alternation of normal files
####################################################
for iter in range(len(names)):
    if update_files[iter] == 0:
        continue
    
    filename = names[iter] + ".html"
    # Read in the file
    with open(filename, 'r') as file:
      filedata = file.read()
    
    # Change header url if necessary
    find_string = "<a class=\"navbar-brand\" href=\"index.html\">"
    replace_string = "<a class=\"navbar-brand\" href=\"" + base_dir[iter] + "index.html\">"
    filedata = filedata.replace(find_string, replace_string)
    
    # Find and elinimate the code for header
    find_string = "<div id=\"header\">\n\n\n\n"
    find_string += "<h1 class=\"title toc-ignore\">"
    find_string += titles[iter]
    find_string += "</h1>\n\n</div>"
    
    replace_string = ""
    filedata = filedata.replace(find_string, replace_string)
    
    # Write the code again
    replace_string = "\n<hr>\n" + find_string + "\n<hr>\n"
    find_string = "<div id=\"content-col\">"
    replace_string = find_string + replace_string
    filedata = filedata.replace(find_string, replace_string)
    
    # Elinimate r codes
    # Initializing a queue
    queue_start = ["1", "2", "3", "4", "5"]
    queue_stop = ["1", "2", "3", "4", "5","6"]
    find_string = "<pre>"
    
    record = False
    stop_loop = False
    # reading each word
    for word in filedata:
        if stop_loop:
            break
        
        if record:
            find_string = find_string + word
            queue_stop.append(word)
            queue_stop.pop(0)
            if queue_stop == ["<", "/", "p", "r", "e", ">"]:
                stop_loop = True                
            
        else:
            queue_start.append(word)
            queue_start.pop(0)
            if queue_start == ["<", "p", "r", "e", ">"]:
                record = True 
    
    replace_string = ""  
    filedata = filedata.replace(find_string, replace_string)                
                

    # Write the file out 
    with open(filename, 'w') as file:
      file.write(filedata)
      
    # Update config file
    config['update'][iter] = 0;

config.to_csv("change_config.csv", index = False)
# ## Read agents result
# setwd("D:/Works/PhD/Research/Processed_Data/csv/output/agent")
# 
# cross_border_data <- read.csv("cross_border.csv", sep = "")
# end_user_data <- read.csv("end_user_type_0.csv", sep = "")
# hydro_HV_data <- read.csv("hydro_HV_plant.csv", sep = "")
# hydro_LV_data <- read.csv("hydro_LV_plant.csv", sep = "")
# wind_HV_data <- read.csv("wind_HV_plant.csv", sep = "")
# wind_LV_data <- read.csv("wind_LV_plant.csv", sep = "")
# 
# ## Read power markets result
# setwd("D:/Works/PhD/Research/Processed_Data/")
# 
# TSO_imbalance_demand <- read.csv("csv/output/power_market/TSO_imbalance_demand.csv", sep = "")

library(readr)
library(rmarkdown)
library("rstudioapi")

# setwd(dirname(getActiveDocumentContext()$path))
# render("index.Rmd", output_format = "html_document")

# setwd(dirname(getActiveDocumentContext()$path))
# render("user_guide.Rmd", output_format = "html_document")

setwd(dirname(getActiveDocumentContext()$path)) 
config <- read.csv("change_config.csv")
render_files <- function(config){
  setwd(dirname(getActiveDocumentContext()$path))
  for(file_iter in 1:nrow(config)){
    if(config$update[file_iter] == 0){
     next
    }
    
    template <- read_lines(paste(config$template[file_iter], ".Rmd", sep = ""))
    content <- read_lines(paste(config$dir[file_iter], config$filename[file_iter], ".Rmd", sep = ""))
    temp <- template
    
    ### Overwrite temp         
    for(line_iter in 1:length(template)){
      ### Change title name
      if(template[line_iter] == "title: \"Template\""){
        temp[line_iter] = paste("title: \"", config$title[file_iter], "\"", sep = "")
      }
      
      ### Change link to current file
      toc_text <- unlist(strsplit(template[line_iter], "[[]"))
      toc_text <- toc_text[length(toc_text)]
      toc_text <- unlist(strsplit(toc_text, "[]]"))
      toc_text <- toc_text[1]
      if(length(toc_text) == 0){
        next
      }
      
      if(toc_text == config$toc_0[file_iter]){
        temp[line_iter] <- paste(unlist(strsplit(template[line_iter], "[(]"))[1], "(#)", sep = "")
        head_tab <- unlist(strsplit(template[line_iter], " "))[1]
      }
      if(toc_text == config$toc_1[file_iter]){
        temp[line_iter] <- paste(unlist(strsplit(template[line_iter], "[(]"))[1], "(#)", sep = "")
      }
      
      ### Insert main content
      if(template[line_iter] == "::: {id=\"content-col\"}"){
        temp[line_iter + 1:length(content)] <- content
        temp <- append(temp, template[(line_iter + 1):length(template)])
        break
      }
    }

    ### create toc
    toc <- c()
    for(line_iter in 1:length(content)){
      if(content[line_iter] == "<!-- title for toc -->"){
        toc_text <- unlist(strsplit(content[line_iter + 1], "[{]"))
        toc_text <- unlist(strsplit(toc_text[1], " "))
        if(length(toc_text) == 2){
          toc_text_temp <- c(paste(head_tab, toc_text[1], sep = ""), paste("[", toc_text[2], "]", sep = ""))
        }else if(length(toc_text) == 3){
          toc_text_temp <- c(paste(head_tab, toc_text[1], sep = ""), paste("[", toc_text[2], sep = ""), paste(toc_text[3], "]", sep = ""))
        }else{
          toc_text_temp <- c(paste(head_tab, toc_text[1], sep = ""), paste("[", toc_text[2], sep = ""), toc_text[3:(length(toc_text) - 1)], paste(toc_text[length(toc_text)], "]", sep = ""))
        }
        toc_text <- toc_text_temp[1]
        for(text_iter in 2:length(toc_text_temp)){
          toc_text <- paste(toc_text, toc_text_temp[text_iter])
        }
        
        toc <- c(toc, toc_text)
      }
    }    
        
    ### Insert intra-page toc 
    for(line_iter in 1:length(temp)){
      toc_text <- unlist(strsplit(temp[line_iter], "[[]"))
      toc_text <- toc_text[length(toc_text)]
      toc_text <- unlist(strsplit(toc_text, "[]]"))
      toc_text <- toc_text[1]
      if(length(toc_text) == 0){
        next
      }      
      
      if(toc_text == config$toc_0[file_iter]){
        temp <- c(temp[1:line_iter], toc, temp[(line_iter + 1):length(temp)])
      }      
    }
    
    write(temp, paste(config$dir[file_iter], "temp.Rmd", sep = ""))
    
    render(paste(config$dir[file_iter], "temp.Rmd", sep = ""), output_file = config$filename[file_iter], output_format = "html_document")
    file.remove(paste(config$dir[file_iter], "temp.Rmd", sep = ""))
  }
}
render_files(config)

# raw_data <- read_lines("examples/template.Rmd")
# write(raw_data, "test.Rmd")
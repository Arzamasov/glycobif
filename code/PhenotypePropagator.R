# load libraries
library(data.table)
library(openxlsx)
library(caret)
library(yaml)

# read command-line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript code/PhenotypePropagator.R <input_list> <output>")
}

magtable_path <- args[1]  # Path to input list file
OUTPUT_DIR <- args[2]   # Output directory

# print arguments for debugging
print(paste("Input List File:", magtable_path))
print(paste("Output Directory:", OUTPUT_DIR))

# set up configuration
config <- yaml.load_file("code/PhenotypePropagator.yaml")

noStyle <- createStyle(fontColour = "#9A0511", bgFill = "#FEC7CE")
yesStyle <- createStyle(fontColour = "#09600B", bgFill = "#C7EECF")

ANNOTATION_DIR <- file.path(config$root_dir,config$annotation_dir)
MODEL_DIR <- file.path(config$root_dir,config$model_dir)
TRAINING_DIR <- file.path(config$root_dir,config$training_dir)
TRAINING_DIR_GENOMEID <- file.path(config$root_dir,config$training_dir_genomeid)

# list of subsystems, phenotypes and training filenames
subsRefPATH <- file.path(config$root_dir,config$subs_ref)
subsRef <- read.csv(subsRefPATH,sep='\t')
subsRef <- data.table(subsRef)
subsRef <- subsRef[isDone == 1]
subsRef[,isDone := NULL]

# list of functional roles
funcRolesPATH <- file.path(config$root_dir,config$functional_roles)
funcRoles <- read.csv(funcRolesPATH,sep='\t')
funcRoles <- data.table(funcRoles)

# list of genomes
refMAGtablePATH <- file.path(magtable_path)
refMAGtable <- read.csv(magtable_path,sep='\t',colClasses = "character")
refMAGtable <- data.table(refMAGtable)

wbMAG <- createWorkbook()
wbFull <- createWorkbook()

results <- data.table()

tableMatML <- data.table()
for(p in 1:nrow(subsRef)){
  curPhenotype <- subsRef[p]$Phenotype
  curSubsystem <- subsRef[p]$Subsystem
  curFilename <- subsRef[p]$Filename

  rolesSubsystem <- funcRoles[subsystem == curSubsystem]
  
  model <- readRDS(file.path(MODEL_DIR,curPhenotype,paste0(curPhenotype,"_ranger.rds")))
  trainset <- read.csv(file.path(TRAINING_DIR,curFilename),colClasses = "factor",check.names = FALSE)
  trainset <- data.table(trainset)
  trainset <- trainset[1]
  trainset[,Y:=NULL]
  for(n in names(trainset)){
    trainset[,(n):="0"]
  }
  
  trainsetGID <- read.csv(file.path(TRAINING_DIR_GENOMEID,curFilename),colClasses = "factor",check.names = FALSE)
  trainsetGID <- data.table(trainsetGID)
 
  files <- list.files(ANNOTATION_DIR)
  
  data4write <- data.table()
  resultsMatML <- c()
  for(f in files){
    print(paste0("Genome: ",f,", Subsystem: ", curSubsystem,", Phenotype: ",curPhenotype))
    
    genome <- gsub(".faa","",f)
        
    data <- read.csv(file.path(ANNOTATION_DIR,f),sep='\t')
    data <- data.table(data)
        
    genePresence <- merge(rolesSubsystem,data,by.x="long_name",by.y="Winner")
    genes <- genePresence$short_name
        
    # prediction by ML ranger
    curset <- copy(trainset)
    for(g in genes){
      if(g %in% names(curset)){
       curset[,(g):="1"]
      }
    }
    set.seed(1)
    phenByMLranger <- predict(model,curset)
    data4write <- rbind(data4write,cbind("ID"=genome,curset))
    
    MAGprops <- refMAGtable[genome_ID == genome]
    MAG_species <- MAGprops[1]$curated_taxonomy
               
    phenByMLranger <- as.integer(as.character(phenByMLranger))
    
    results <- rbind(results,data.table("genome"=genome,
                                        "curated_taxonomy"=MAG_species,
                                        "subsystem"=curSubsystem,
                                        "phenotype"=curPhenotype,
                                        "phenByMLranger"=phenByMLranger
    ))
    
       resultsMatML <- c(resultsMatML,phenByMLranger)
  }
  
  # full output
  if(nrow(results) > 0){ 
   subsResults <- results[phenotype == curPhenotype]
   subsResults[,subsystem:=NULL]
   subsResults[,phenotype:=NULL]
   addWorksheet(wbFull,curPhenotype)
   data4write[,ID:=NULL]
   writeData(wbFull,sheet=p,cbind(subsResults,data4write),startRow=1)
   conditionalFormatting(wbFull,sheet=p,rows=2:(nrow(subsResults)+1),cols=(ncol(subsResults)+1):(ncol(subsResults)+ncol(data4write)),type="contains",rule="1",yesStyle)
   conditionalFormatting(wbFull,sheet=p,rows=2:(nrow(subsResults)+1),cols=(ncol(subsResults)+1):(ncol(subsResults)+ncol(data4write)),type="contains",rule="0",noStyle)
  }
  
  # BPM
  if(nrow(results) > 0){
   tableMatML <- cbind(tableMatML,"tmp"=resultsMatML)
   setnames(tableMatML,"tmp",curPhenotype)
  }
}

saveWorkbook(wbFull,file.path(OUTPUT_DIR,"PredictionsFull.xlsx"), overwrite = TRUE)

tableMatML <- cbind("genome"=subsResults$genome,
                    "curated_taxonomy"=subsResults$curated_taxonomy,
                     tableMatML)

wbMat <- createWorkbook()
addWorksheet(wbMat,"Predictions")
writeData(wbMat,sheet=1,tableMatML,startRow=1)
saveWorkbook(wbMat,file.path(OUTPUT_DIR,"BPM.xlsx"), overwrite = TRUE)

# 
# ############################################################
# #                                                          #
# ###########################HUMAN############################
# #                                                          #
# ############################################################
# 
# 
# # By this script I want to creat known structure PPI dataset of MIto proteome  
# 
# library(vanddraabe)
# library(rgl)
# library(Rpdb)
library(bio3d)
# library(reshape2)
library(ggplot2)
# library(seqinr)
# require("Biostrings")
# #data(BLOSUM100)
# library(msa)
# library(RColorBrewer)
# library(gridExtra)
# library(ggpubr)
 library(e1071)     
# library(fitdistrplus)
# library(pryr)
# library(pracma) 
 #library(dplyr) #for the revalue function  
 library(tidyr)
# library(ggjoy)
# library(plotrix)
# library(ggExtra)
library(pracma)
#library(reshape2)
working_dir="/media/elhabashy/Elements/Hadeer_backup/hybrid_Xlmsecs"

#human
list0=data.frame(read.delim(paste(working_dir,"about_dataset","mito_proteome/published_human_mitocarta3.0_impiq2_proteins.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
list=list0

#Creata blast file for all proteins
for(i in 1:nrow(list)){
  #for(i in 1:10){
  tryCatch({
    print(i)
    
    #create a directory for the every protein
    protein_dir=" "
    protein_dir=paste(working_dir,"dataset/mito_proteome/human",list[i,"uid"],sep="/")
    dir.create(protein_dir)
    
    
    #download protein info from the UniProt
    #download.file(paste("https://www.uniprot.org/uniprot/?query=",list[i,"uid"],"&sort=score&columns=id,entry%20name,protein%20names,genes,genes(ALTERNATIVE),organism,length,3d,interactor,comment(ALTERNATIVE%20PRODUCTS),database(Pfam)&format=tab",sep=""), destfile = paste(protein_dir, paste(list[i,"uid"],"_info.csv",sep = ""), sep='/'))
    
    #download protein fasta sequences from the UniProt
    #download.file(paste("https://www.uniprot.org/uniprot/",list[i,"uid"],".fasta",sep=""), destfile = paste(protein_dir, paste(list[i,"uid"],".fasta",sep = ""), sep='/'))
    
    
    # this has been done to 
    #blastp –query file.fasta –db pdbaa –out output.txt -evalue
    #system(paste("cd",  protein_dir, "&& /Hadeer/software/ncbi-blast-2.9.0+/bin/blastp -query",  paste(list[i,"uid"],".fasta",sep = "") , "-db /Hadeer/software/ncbi-blast-2.9.0+/blastdb/pdbaa/pdbaa","-evalue 0.00001", "-out" ,paste(list[i,"uid"], "blastp_feb2023.out",sep="_"), "-outfmt 6", sep=" "))
    system(paste("cd",  protein_dir, "&& /Hadeer/software/ncbi-blast-2.9.0+/bin/blastp -query",  paste(list[i,"uid"],".fasta",sep = "") , "-db /Hadeer/software/ncbi-blast-2.9.0+/blastdb/pdbaa_28_10_2024/pdbaa","-evalue 0.00001", "-out" ,paste(list[i,"uid"], "blastp_oct2024.out",sep="_"), "-outfmt 6", sep=" "))
    
    #-outfmt " Score evalue Identities Gaps "
    # /Hadeer/software/ncbi-blast-2.9.0+/blastdb/pdbaa
    #dataset 
    
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}



struc_ppi=data.frame()
#Creata blast file for all proteins
for(i in 1:nrow(list)){
  #for(i in 1:20){
  tryCatch({
    print(i)
    
    #create a directory for the every protein
    protein1_dir=" "
    protein1_dir=paste(working_dir,"dataset/mito_proteome/human",list[i,"uid"],sep="/")
    #dir.create(protein_dir)
    
    temp_file1=data.frame()
    temp_file1=data.frame(read.delim(paste(protein1_dir,paste(list[i,"uid"], "blastp_oct2024.out",sep="_"), sep='/'), header=FALSE, stringsAsFactors = FALSE))
    names(temp_file1)[1:12]= c("query", "subject", "%identity", "alignment length", "mismatches", "gap opens", "q. start:", "q. end", "s. start", "s. end", "evalue", "bit score")
    temp_file1[,"pdb"]=gsub("_.*", "\\4", temp_file1[,"subject"])
    temp_file1[,"chain"]=gsub(".*_", "\\1", temp_file1[,"subject"])
    temp_file1=temp_file1[!duplicated(temp_file1[,"subject"]),]
    temp_file1=temp_file1[temp_file1[,"%identity"]>=50,]
    
    for(j in 1:nrow(list)){
      tryCatch({
        print(paste(i,"=>",j))
        
        #create a directory for the every protein
        protein2_dir=" "
        protein2_dir=paste(working_dir,"dataset/mito_proteome/human",list[j,"uid"],sep="/")
        #dir.create(protein_dir)
        
        temp_file2=data.frame()
        temp_file2=data.frame(read.delim(paste(protein2_dir,paste(list[j,"uid"], "blastp_oct2024.out",sep="_"), sep='/'), header=FALSE, stringsAsFactors = FALSE))
        names(temp_file2)[1:12]= c("query", "subject", "%identity", "alignment length", "mismatches", "gap opens", "q. start:", "q. end", "s. start", "s. end", "evalue", "bit score")
        temp_file2[,"pdb"]=gsub("_.*", "\\4", temp_file2[,"subject"])
        temp_file2[,"chain"]=gsub(".*_", "\\1", temp_file2[,"subject"])
        temp_file2=temp_file2[!duplicated(temp_file2[,"subject"]),]
        temp_file2=temp_file2[temp_file2[,"%identity"]>=50,]
        
        #remove case where 
        temp_file=data.frame()
        temp_file=merge(temp_file1,temp_file2,by.x="pdb", by.y="pdb")
        temp_file= temp_file[temp_file[,"subject.x"]!=temp_file[,"subject.y"],]   
        
        #remove cases where 
        temp_file[,"test1"]=paste(temp_file[,"subject.x"], temp_file[,"subject.y"],sep=" ")
        temp_file[,"test2"]=paste(temp_file[,"subject.y"], temp_file[,"subject.x"],sep=" ")
        
        temp_file[,"delete"]=" "
        for(k in 1:nrow(temp_file)){
          for(l in 1:nrow(temp_file)){
            if(temp_file[k,"test1"]==temp_file[l,"test2"]){
              temp_file[k,"delete"]=TRUE
              temp_file[l,"delete"]=TRUE
              print("case detected")
            }
          }
        } 
        
        temp_file= temp_file[temp_file[,"delete"]!=TRUE,]   
        
        
        if(nrow(temp_file)!=0){
          
          write.csv( temp_file, file =paste("/Hadeer/hybrid_Xlmsecs/dataset/mito_complex/human/",list[i,"uid"],"-",list[j,"uid"],"_oct2024.csv",sep=""))

          struc_ppi_temp=data.frame()
          struc_ppi_temp[1,1:5]=list[i,c("Entry","uid","Organism","Length","Cross.reference..Pfam.")]
          names(struc_ppi_temp)[1:5]=c("interactorA_entry","interactorA_uid", "interactorA_organism","interactorA_length","interactorA_Pfam_domains")
          
          struc_ppi_temp[1,6:10]=list[j,c("Entry","uid","Organism","Length","Cross.reference..Pfam.")]
          names(struc_ppi_temp)[6:10]=c("interactorB_entry","interactorB_uid", "interactorB_organism","interactorB_length","interactorB_Pfam_domains")
          struc_ppi_temp[,"complex_3d"]=as.numeric(nrow(temp_file))
          
          struc_ppi=rbind(struc_ppi,struc_ppi_temp)
        }
        
      }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
      
    }
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}



##############################################################
##############################################################
##############################################################
#Add prefixes 
list_backup=list
list=struc_ppi

#Add prefixes 
#list0=data.frame(read.delim("/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_positive_mito_ppi.csv",header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
#list=list0
list[,"uid1"]=list[,"interactorA_uid"]
list[,"uid2"]=list[,"interactorB_uid"]
list[,"test"]=paste(list[,"interactorA_uid"],list[,"interactorB_uid"],sep=" ")
list[,"test2"]=paste(list[,"interactorB_uid"],list[,"interactorA_uid"],sep=" ")
#/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mitocomplex_p6_prefix.csv
#list of alignable pairs 
list1=data.frame(read.delim(paste(working_dir,"about_dataset/mito_proteome/human_mitocomplex_p1_prefix.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
list2=data.frame(read.delim(paste(working_dir,"about_dataset/mito_proteome/human_mitocomplex_p2_prefix.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
list3=data.frame(read.delim(paste(working_dir,"about_dataset/mito_proteome/human_mitocomplex_p3_prefix.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
list4=data.frame(read.delim(paste(working_dir,"about_dataset/mito_proteome/human_mitocomplex_p4_prefix.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
list5=data.frame(read.delim(paste(working_dir,"about_dataset/mito_proteome/human_mitocomplex_p5_prefix.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
list6=data.frame(read.delim(paste(working_dir,"about_dataset/mito_proteome/human_mitocomplex_p6_prefix.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))

list_all=rbind(list1,list2,list3,list4,list5,list6)
rm(list1,list2,list3,list4,list5,list6)
list_all[,"test"]=paste(list_all[,"uid1"],list_all[,"uid2"],sep=" ")


for(i in 1:nrow(list)){
  print(i)
  list_all_temp=data.frame()
  list_all_temp=list_all[list_all[,"test"]==list[i,"test"],]
  if(nrow(list_all_temp)==1){ 
    list[i,16:40]=  list_all_temp[1,1:25]
  }
}


for(i in 1:nrow(list)){
  print(i)
  
  list_all_temp=data.frame()
  list_all_temp=list_all[list_all[,"test"]==list[i,"test2"],]
  if(nrow(list_all_temp)==1){ 
    list[i,16:40]=  list_all_temp[1,1:25]
  }
}


list=list[,-c(1:15)]
list=list[!is.na(list[,"prefix"]),]
list=list[!duplicated(list[,"prefix"]),]
write.csv(list, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d_6245_Nov2024.csv", quote = TRUE)
names(list)[4]="uid1"
names(list)[14]="uid2"


#########################################
#calculate pair_wise seqidentity  
#########################################
#global
#Run pair wise seqeunce similarity between pairs
for(i in 1:nrow(list)){
  tryCatch({
    print(i)
    #if(is.na(list[i,"complex_pdb"])){
    #create a directory for the every protein
    protein1_dir=" "
    protein1_dir=paste(working_dir,"dataset/mito_proteome/human",list[i,"uid1"],sep="/")
    
    #create a directory for the every protein
    protein2_dir=" "
    protein2_dir=paste(working_dir,"dataset/mito_proteome/human",list[i,"uid2"],sep="/")
    
    #create a directory for the every protein
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset/mito_complex/human/complexes", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
    dir.create(complex_dir)
    
    #dir.create(complex_dir)   
    #temp_dir=" "
    #temp_dir=paste("/media/hadeer/Elements/human_mito/benckmark/positive/",sep = "")
    
    # this has been done to 
    #blastp –query file.fasta –db pdbaa –out output.txt -evalue
    #system(paste("cd",temp_dir, "&& /Hadeer/software/emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0/emboss/needle",  "-asequence", paste(protein1_dir, paste(list[i,"uid1"],".fasta",sep = ""), sep='/') , "-bsequence", paste(protein2_dir, paste(list[i,"uid2"],".fasta",sep = ""), sep='/'), "-out" ,paste(list[i,"uid1"],"vs",list[i,"uid2"], "needle.out",sep="_"), sep=" "))
    system(paste("cd", complex_dir, "&& /Hadeer/software/emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0/emboss/stretcher",  "-asequence", paste(protein1_dir, paste(list[i,"uid1"],".fasta",sep = ""), sep='/') , "-bsequence", paste(protein2_dir, paste(list[i,"uid2"],".fasta",sep = ""), sep='/'), "-out" ,paste(list[i,"uid1"],"vs",list[i,"uid2"], "stretcher.out",sep="_"), sep=" "))
    
    temp_file=data.frame()
    temp_file=data.frame(read.delim(paste(complex_dir,paste(list[i,"uid1"],"vs",list[i,"uid2"], "stretcher.out",sep="_"), sep='/'), header=FALSE, stringsAsFactors = FALSE))
    temp_file=data.frame(grep("# Identity: ", temp_file[,"V1"],value = TRUE))
    list[i,"pairwise_identity"]=gsub(".*# Identity:    ", "", temp_file[,1])
    system(paste("cd",complex_dir, "&& rm " ,paste(list[i,"uid1"],"vs",list[i,"uid2"], "stretcher.out",sep="_"), sep=" "))
    system(paste(" rm -r " ,complex_dir, sep=" "))
  #  }   
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
} 


list[,"pairwise_identity"]=gsub(".*\\(", "", list[,"pairwise_identity"])
list[,"pairwise_identity"]=gsub("\\).*", "", list[,"pairwise_identity"])
list[,"pairwise_identity"]=gsub("\\%.*", "", list[,"pairwise_identity"])



###################################################

#dirty interface calculation
# I will try to check if there is interface in the first pdb found by blast
for(i in 1:nrow(list)){
  #for(i in 1:20){
  tryCatch({
    #print(i)
   # if(is.na(list[i,"complex_pdb"])){
      print(i)
      #create a directory for the every protein
      complex_dir=" "
      complex_dir=paste("/media/hadeer/Elements/human_mito/benckmark/positive/",paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep="_"),sep="")
      dir.create(complex_dir)
      
      complex_pdb=" "
      complex_pdb=paste(complex_dir,"complex_pdb",sep="/")
      dir.create(complex_pdb)
      
      #some times pdb has many copies of the interactions
      temp_file=data.frame()
      temp_file=data.frame(read.delim(paste("/Hadeer/hybrid_Xlmsecs/dataset/mito_complex/human",paste(list[i,"uid1"],"-",list[i,"uid2"],"_oct2024.csv",sep=""),sep="/"), header=TRUE,sep=",",as.is=TRUE, colClasses="character", stringsAsFactors = FALSE))
      temp_file=temp_file[temp_file[,"pdb"]==temp_file[1,"pdb"],]
      temp_file[,"no_pairs"]=as.numeric(nrow(temp_file))
      download.file(paste(paste("https://files.rcsb.org/view", temp_file[1,"pdb"], sep = '/'), "cif", sep = '.'), destfile = paste(complex_pdb,"/", temp_file[1,"pdb"],".cif",sep=""))
      #get.pdb(temp_file[1,"pdb"], path =paste(complex_pdb,"/", temp_file[1,"pdb"],".cif",sep=""), URLonly=T)
      pdb_file=data.frame()
      pdb_file= read.cif(paste(complex_pdb,"/", temp_file[1,"pdb"],".cif",sep=""))
      
      #pdb_file= read.pdb(paste(complex_pdb,"/", temp_file[1,"pdb"],".pdb",sep=""))
      if(nrow(pdb_file$atom)==0){next}
      chain1 =data.frame()
      chain1 = trim.pdb( pdb_file, chain=as.character(temp_file[1,"chain.x"]),type = 'ATOM' )
      #write.pdb(chain1, type = 'ATOM', file = paste(complex_pdb, paste(temp_file[1,"pdb"],"_",temp_file[1,"chain.x"],".pdb",sep = ""), sep = '/'))
      
      chain2 =data.frame()
      chain2 = trim.pdb( pdb_file, chain=as.character(temp_file[1,"chain.y"]),type = 'ATOM')
      #write.pdb(chain2,type = 'ATOM',  file = paste(complex_pdb, paste(temp_file[1,"pdb"],"_",temp_file[1,"chain.y"],".pdb",sep = ""), sep = '/'))
      
      chain=cat.pdb( chain1,  chain2, rechain=FALSE)
      #write.pdb(chain,type = 'ATOM',  file = paste(complex_pdb, paste(temp_file[1,"pdb"],"_",temp_file[1,"chain.x"],"_",temp_file[1,"chain.y"],".pdb",sep = ""), sep = '/'))
      
      #some times
      
      for(l in 1:nrow(temp_file)){
        na_interface=data.frame()
        chain1=data.frame()
        chain1 = trim.pdb( pdb_file, chain=as.character(temp_file[l,"chain.x"]) ,type = 'ATOM')
        chain1 = trim.pdb( chain1, elety="CA")
        
        #chain1=chain1[as.character(chain1$atom[,"elety"])=="CA",]
        chain2=data.frame()
        chain2 = trim.pdb(pdb_file, chain=as.character(temp_file[l,"chain.y"]) ,type = 'ATOM')
        chain2 = trim.pdb( chain2, elety="CA")
        #chain2=atom.select(chain2, "calpha")
        
        
        for (j in 1:nrow(chain1$atom)) {
          #for (j in 1:30) {
          Xj= chain1$atom[j,"x"]
          Yj= chain1$atom[j,"y"]
          Zj= chain1$atom[j,"z"]
          
          for (k in 1:nrow(chain2$atom)){  
            Xk= chain2$atom[k,"x"]
            Yk= chain2$atom[k,"y"]
            Zk= chain2$atom[k,"z"]
            
            dista= round( sqrt((Xj-Xk)^2+(Yj-Yk)^2+(Zj-Zk)^2),digits = 4)
            
            if(dista <= 10){
              na_interface_temp=data.frame(test=c(" "))
              na_interface_temp[,"test"]=paste( chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],sep=" ")
              na_interface_temp[,"ele"]=paste(chain1$atom[j,"elety"],chain2$atom[k,"elety"], sep=" ")
              na_interface_temp[,"na_dist"]=paste(dista)
              na_interface=rbind(na_interface,na_interface_temp)}
          }
          
          
          if(nrow(na_interface)>10){break}
        }
        
        if(nrow(na_interface)>0){ temp_file[l,"interface"]=TRUE}else{temp_file[l,"interface"]=FALSE}
        
        #complex_dir=" "
        #complex_dir=paste(working_dir,"dataset/complexes", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
        #dir.create(complex_dir)
        
        #save the complex pair
        # chain1=data.frame()
        # chain1 = trim.pdb( pdb_file, chain=as.character(temp_file[l,"chain.x"]) ,type = 'ATOM')
        # chain2=data.frame()
        # chain2 = trim.pdb(pdb_file, chain=as.character(temp_file[l,"chain.y"]) ,type = 'ATOM')
        # chain12=cat.pdb(chain1, chain2, rechain=FALSE)
        # chain12= clean.pdb(chain12, consecutive =FALSE, force.renumber = FALSE, fix.chain= FALSE, fix.aa= TRUE, rm.wat= TRUE, rm.lig= TRUE, rm.h= TRUE, verbose= TRUE)
        # complex_pdb=" "
        # complex_pdb=paste(complex_dir, "complex_pdbs",sep="/")
        # dir.create(complex_pdb)
        # write.pdb(chain12,type = 'ATOM', file= paste(complex_pdb, paste(paste(temp_file[l,"pdb"],temp_file[l,"chain.x"],temp_file[l,"chain.y"],sep='_'), "pdb", sep = '.'), sep='/'))
        
      }
      
      
      list[i,"interface"]=as.numeric(nrow(temp_file[temp_file[,"interface"]==TRUE,]))
      
      temp_file1=data.frame          
      temp_file1=temp_file[temp_file[,"interface"]==TRUE & !duplicated(temp_file[,"chain.x"]),]
      temp_file2=data.frame 
      temp_file2=temp_file[temp_file[,"interface"]==TRUE & !duplicated(temp_file[,"chain.y"]),]  
      temp_file[1,"interface_architect"]=paste(nrow(temp_file1), ":" ,nrow(temp_file2), sep=" ")
      
      
      list[i,"complex_pdb"]=temp_file[1,"pdb"]
      list[i,"complex_chain1"]=temp_file[1,"chain.x"]
      list[i,"complex_chain2"]=temp_file[1,"chain.y"]
      list[i,"no_pairs"]=temp_file[1,"no_pairs"]
      list[i,"interface_architect"]=temp_file[1,"interface_architect"]
      
      system(paste("rm ", paste(complex_pdb,"/", temp_file[1,"pdb"],".cif",sep="") ))
   # }
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}
write.csv(list, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d_6245_Nov2024.csv", quote = TRUE)
list_all=list
list_interface=list[list[,"interface"]==1,]
list_interface=list_interface[!duplicated(list_interface[,"prefix"]),]
list_interface=list_interface[!is.na(list_interface[,"ecs_complete"]),]
list_interface=list_interface[!is.na(list_interface[,"kurtosis"]),]
#################################
#check for complete ecs
list=data.frame(read.delim(paste(working_dir,"about_dataset","mito_proteome/human_mito_structure_ppi_prefix_3d_6245_Nov2024.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))

set=data.frame(v1=c("00","01","02","03","04","05","06","07","08","09","10","11","12","13","14"))
list[,"ecs_complete"]= " "
list[,"ecs_anomaly"]= " "

for (i in 1:nrow(list)){
  #for (i in 1:50){
  for (k in 1:nrow(set)){
    tryCatch({
      
      print(i)
      ec_list=data.frame()
      ec_list=data.frame(read.delim(paste("/media/hadeer/Elements/human_mito/evcomplex",paste("sepmito_p",set[k,"v1"],sep=""), paste(list[i,"prefix"], "CouplingScores_inter.csv", sep='_'), sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
      if(nrow(ec_list)!=0){list[i,"ecs_complete"]="Complete"}
      if(ec_list[1,"A_i"]=="-"){list[i,"ecs_anomaly"]="TRUE"}
      
    },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
  }
}



for (i in 1:nrow(list)){
  for (k in 1:nrow(set)){
    tryCatch({
    
    print(i)
    ec_list=data.frame()
    ec_list=data.frame(read.delim(paste("/media/hadeer/Elements/human_mito/evcomplex",paste("sepmito_p",set[k,"v1"],sep=""), paste(list[i,"prefix"], "CouplingScores_inter.csv", sep='_'), sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    
    if(nrow(ec_list)!=0){list[i,"ecs_complete"]="Complete"
    
    ec_list= ec_list[order( as.numeric(ec_list[,"cn"]), decreasing = TRUE),]
    
    list[i,"kurtosis"]=kurtosis(as.numeric(ec_list[,"cn"]))
    list[i,"skewness"]=skewness(as.numeric(ec_list[,"cn"]))
    list[i,"squared_skewness"]=( list[i,"skewness"])^2
    list[i,"max_ecs"]=max(as.numeric(ec_list[,"cn"]))
    }
    
    
   },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
  }
}

#########################################

list[,"af3"]=NA
list[,"test"]=paste(list[,"prefix"],sub("\\_.*", "",list[,"uid1"]), sub("\\_.*", "",list[,"uid2"]),sep='_')


for(i in 1:nrow(list)){
  tryCatch({
    print(i)
    af3_list=data.frame()
    af3_list=data.frame(read.delim(paste("/media/hadeer/Elements/human_mito/alphafold3/af3_list.txt"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    af3_list=data.frame(af3_list[af3_list[,1]==tolower(list[i,"test"]),])
    
    if(nrow(af3_list)!=0){list[i,"af3"]="complete"
    
    af_score=data.frame()
    af_score=jsonlite::fromJSON(readLines(paste("/media/hadeer/Elements/human_mito/alphafold3",paste("fold",list[i,"prefix"],sub("\\_.*", "",tolower(list[i,"uid1"])), sub("\\_.*", "",tolower(list[i,"uid2"])),sep='_'),paste("fold",list[i,"prefix"],sub("\\_.*", "",tolower(list[i,"uid1"])), sub("\\_.*", "",tolower(list[i,"uid2"])),"summary_confidences_0.json",sep='_'),sep="/")),flatten=TRUE)
    af_score_temp=as.data.frame( af_score)
    list[i,"iptm"]=af_score_temp[1,"iptm"]
    list[i,"ptm"]=af_score_temp[1,"ptm"]
    list[i,"ranking_score"]=af_score_temp[1,"ranking_score"]
    
    } 
    
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
}

list_all=list
list_interface=list[list[,"interface"]==1,]
list_interface=list_interface[!duplicated(list_interface[,"prefix"]),]
list_interface=list_interface[!is.na(list_interface[,"ecs_complete"]),]
list_interface=list_interface[!is.na(list_interface[,"kurtosis"]),]
list=list_interface
write.csv(list, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d_714_Nov2024.csv", quote = TRUE)

list_af3=list[is.na(list[,"af3"]),]
write.csv(list_af3, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d_714_Nov2024_torunaf3.csv", quote = TRUE)

#################################
#check diagonal beahvior
#plot contacts of ECs


for (i in 1:nrow(list)){
  for (k in 1:nrow(set)){
    tryCatch({
    print(i)
    ec_list=data.frame()
    ec_list=data.frame(read.delim(paste("/media/hadeer/Elements/human_mito/evcomplex",paste("sepmito_p",set[k,"v1"],sep=""), paste(list[i,"prefix"], "CouplingScores_inter.csv", sep='_'), sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    
    if(nrow(ec_list)!=0){list[i,"ecs_complete"]="Complete"
    ec_list=ec_list[order(ec_list[,"cn"],decreasing=T),]
    ec_list=head(ec_list,50)
    
    EC_contact=ggplot(ec_list, aes(x=ec_list[,"i"],y=ec_list[,"j"],color=ec_list[,"cn"]))
    EC_contact=EC_contact+ geom_point()
    EC_contact=EC_contact+ xlim(1,max(ec_list[,"i"]))+ ylim(1,max(ec_list[,"j"]))
    EC_contact=EC_contact+  xlab("Residues of protein 1")+theme_bw()
    EC_contact=EC_contact+ ylab("Residues of protein 2")
    EC_contact
    ggsave(paste("/media/hadeer/Elements/human_mito/evcomplex/positive/figures/diagonal", paste(list[i,"prefix"], "CouplingScores_inter_string.png", sep='_'), sep="/"),  scale = 1, width = NA, height = NA, dpi=300)
    }
    },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
  }
}



list[,"diagonal"]=NA

for(i in 1:nrow(list)){
  tryCatch({
    print(i)
    
    diag_list=data.frame()
    diag_list=data.frame(read.delim(paste("/media/hadeer/Elements/human_mito/evcomplex/positive/figures/diagonal/diag.txt",sep="/"),header=F, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    diag_list=data.frame(diag_list[diag_list[,"V1"]==list[i,"prefix"],])
    
    if(nrow(diag_list)!=0){
      print("case detected")
      list[i,"diagonal"]="TRUE"}else{list[i,"diagonal"]="FALSE"}
    
    
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}
###################################


#check which one of those are known in string 
#with conditions of score 0.4 and physical interaction

for (i in 1:nrow(list)){
  # for (i in 1:10){
  
  tryCatch({
    
    #create a directory for the every protein
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset/mito_complex/human/complexes", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
    dir.create(complex_dir)
    
    download.file(paste("https://string-db.org/api/tsv/network?identifiers=", list[i,"uid1"] ,"%0d",list[i,"uid2"],"&required_score=400&network_type=physical",sep=""), destfile = paste(complex_dir, paste(list[i,"uid1"] ,"-",list[i,"uid2"],"_string.csv",sep = ""), sep='/'))
    string=data.frame()
    string=data.frame(read.delim( paste(complex_dir, paste(list[i,"uid1"] ,"-",list[i,"uid2"],"_string.csv",sep = ""), sep='/'),header=TRUE, row.name=NULL, sep="\t",comment.char='#',stringsAsFactors = FALSE))
    list[i,"string"]=nrow(string)
    system(paste(" rm -r " ,complex_dir, sep=" "))
    
  },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}

####################################
###check linearity


set=data.frame(v1=c("00","01","02","03","04","05","06","07","08","09","10","11","12","13","14"))
list[,"R_squared"]= " "
list[,"r"]=" "
list[,"R2_predict_diag"]=" "
list[,"r_predict_diag"]=" "

for (i in 1:nrow(list)){
  #for (i in 1:50){
  for (k in 1:nrow(set)){
    tryCatch({
      
      print(i)
      ec_list=data.frame()
      ec_list=data.frame(read.delim(paste("/media/hadeer/Elements/human_mito/evcomplex",paste("sepmito_p",set[k,"v1"],sep=""), paste(list[i,"prefix"], "CouplingScores_inter.csv", sep='_'), sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
      if(nrow(ec_list)!=0){list[i,"ecs_complete"]="Complete"
      ec_list=ec_list[order(ec_list[,"cn"],decreasing=T),]
      ec_list=head(ec_list,100)

      #2. R SQUARED error metric -- Coefficient of Determination
      ml = lm(i~j, data = ec_list) 
      
      # Extracting R-squared parameter from summary 
      list[i,"R_squared"]=summary(ml)$r.squared
      list[i,"r"]=cor(ec_list[,"i"], ec_list[,"j"], method = 'pearson')
      if(as.numeric(list[i,"R_squared"])>=0.9){list[i,"R2_predict_diag"]="TRUE"}
      if(as.numeric(list[i,"R_squared"])<0.9){list[i,"R2_predict_diag"]="FALSE"}
      
      if(as.numeric(list[i,"r"])>=0.9){list[i,"r_predict_diag"]="TRUE"}
      if(as.numeric(list[i,"r"])<0.9){list[i,"r_predict_diag"]="FALSE"}
      }
      
      
    },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
  }
}

conf_matrix=data.frame(table(True=list[,"diagonal"],Predicted=list[,"R2_predict_diag"]))
cm_plot= ggplot(conf_matrix, aes(x=conf_matrix[,"True"] , y = conf_matrix[,"Predicted"] , fill=conf_matrix[,"Freq"])  ) 
cm_plot= cm_plot+ scale_fill_gradientn( name = "Count", colors = c("white","#B4C79C","darkolivegreen3", "mediumpurple3"))
  #scale_fill_gradient(aesthetics = "fill",name = "Count", palette = pal_seq_gradient(low = c("white", "#56941e", "#7b52ae")))
cm_plot= cm_plot+ geom_text(aes(label = conf_matrix[,"Freq"]), color = "black", size = 6, fontface = "bold")
#scale_colour_gradient( name = waiver(),low = "#94c773",high = "#7b52ae",space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour" )
cm_plot= cm_plot+ geom_tile(color = "white", lwd = 1.5, linetype = 1) +  coord_fixed()
cm_plot= cm_plot+labs( title = "Confusion Matrix", x = "Actual", y = "Predicted") 
cm_plot= cm_plot+geom_text(aes(label = conf_matrix[,"Freq"]), color = "black", size = 6, fontface = "bold", vjust = 0.5, hjust = 0.5) 
cm_plot=cm_plot+ theme_bw()+theme(plot.title = element_text(hjust = 0.5, size = 16),axis.text = element_text(size = 12),  axis.title = element_text(size = 14))
cm_plot




conf_matrix=data.frame(table(True=list[,"diagonal"],Predicted=list[,"r_predict_diag"]))
cm_plot= ggplot(conf_matrix, aes(x=conf_matrix[,"True"] , y = conf_matrix[,"Predicted"] , fill=conf_matrix[,"Freq"])  ) 
cm_plot= cm_plot+ scale_fill_gradientn( name = "Count", colors = c("white","#B4C79C","darkolivegreen3", "mediumpurple3"))
#scale_fill_gradient(aesthetics = "fill",name = "Count", palette = pal_seq_gradient(low = c("white", "#56941e", "#7b52ae")))
cm_plot= cm_plot+ geom_text(aes(label = conf_matrix[,"Freq"]), color = "black", size = 6, fontface = "bold")
#scale_colour_gradient( name = waiver(),low = "#94c773",high = "#7b52ae",space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour" )
cm_plot= cm_plot+ geom_tile(color = "white", lwd = 1.5, linetype = 1) +  coord_fixed()
cm_plot= cm_plot+labs( title = "Confusion Matrix", x = "Actual", y = "Predicted") 
cm_plot= cm_plot+geom_text(aes(label = conf_matrix[,"Freq"]), color = "black", size = 6, fontface = "bold", vjust = 0.5, hjust = 0.5) 
cm_plot=cm_plot+ theme_bw()+theme(plot.title = element_text(hjust = 0.5, size = 16),axis.text = element_text(size = 12),  axis.title = element_text(size = 14))
cm_plot





ggplot(list, aes(x= list[,"r"], y=list[,"R_squared"],color=list[,"diagonal"])) + geom_point()
library(dplyr)
colors = c("white","#B4C79C", "#56941e","#B09FCA", "#7b52ae"))
conf_matrix <- list %>%
  mutate(
    true_label = ifelse(diagonal, "TRUE", "FALSE"),
    pred_label = ifelse(predict_diag, "TRUE", "FALSE")
  ) %>%
  count(true_label, pred_label) %>%
  complete(true_label, pred_label, fill = list(n = 0))

# Step 4: Visualize confusion matrix with ggplot2
ggplot(conf_matrix, aes(x = true_label, y = pred_label, fill = n)) +
  geom_tile(color = "black") +
  geom_text(aes(label = n), color = "white", size = 6) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(
    title = "Confusion Matrix",
    x = "Actual",
    y = "Predicted",
    fill = "Count"
  ) +
  theme_minimal()
rsquaredplot= ggplot(list, aes(x= list[,"r"], y=list[,"R_squared"],color=list[,"diagonal"])) + geom_point()

diag_threshould=0.9
# Step 2: Predict classes based on R_squared


list[,"predicted"]=list[list[,"R_squared"]>0.9,]
library(dplyr)

conf_matrix <- data %>%
  mutate(
    true_label = ifelse(diagonal, "True", "False"),
    pred_label = ifelse(predicted, "True", "False")
  ) %>%
  count(true_label, pred_label) %>%
  complete(true_label, pred_label, fill = list(n = 0))


EC_contact=ggplot(ec_list, aes(x=ec_list[,"i"],y=ec_list[,"j"],color=ec_list[,"cn"]))
EC_contact=EC_contact+ geom_point()
EC_contact=EC_contact+ xlim(1,max(ec_list[,"i"]))+ ylim(1,max(ec_list[,"j"]))
EC_contact=EC_contact+  xlab("Residues No. of protein 1")+theme_bw()
EC_contact=EC_contact+ ylab("Residues No. of protein 2")
EC_contact

write.csv(list, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d_714_Nov2024.csv", quote = TRUE)

###################################
#comparing strings of Pfam domains to find common domains

for(i in 1:nrow(list)){
  #for(i in 1:10){
  tryCatch({
    print(i)
    #create a directory for the every protein
    protein1_dir=" "
    protein1_dir=paste(working_dir,"dataset/mito_proteome/human",list[i,"uid1"],sep="/")
    #dir.create(protein_dir)
    
    #download protein info from the UniProt
    download.file(paste("https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cgene_names%2Clength%2Ccc_subcellular_location%2Cxref_interpro%2Cstructure_3d&format=tsv&query=",list[i,"uid1"],"&format=tab",sep=""), destfile = paste(protein1_dir, paste(list[i,"uid1"],"_info.csv",sep = ""), sep='/'))
    
    temp_file1=data.frame()
    temp_file1=data.frame(read.delim(paste(protein1_dir, paste(list[i,"uid1"],"_info.csv",sep = ""),sep="/"), sep="\t",header=T, stringsAsFactors = FALSE))
    temp_file1=temp_file1[temp_file1[,"Entry.Name"]==list[i,"uid1"],]
    list[i,"uid1_Pfam"]=temp_file1[1,"InterPro"]
    list[i,"uid1_scl"]=temp_file1[1,5]
    
    system(paste("rm ",  paste(protein1_dir, paste(list[i,"uid1"],"_info.csv",sep = ""), sep='/')))
    
    
    
    #create a directory for the every protein
    protein2_dir=" "
    protein2_dir=paste(working_dir,"dataset/mito_proteome/human",list[i,"uid2"],sep="/")
    #dir.create(protein_dir)
    
    #download protein info from the UniProt
    download.file(paste("https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cgene_names%2Clength%2Ccc_subcellular_location%2Cxref_interpro%2Cstructure_3d&format=tsv&query=",list[i,"uid2"],"&format=tab",sep=""), destfile = paste(protein2_dir, paste(list[i,"uid2"],"_info.csv",sep = ""), sep='/'))
    
    temp_file2=data.frame()
    temp_file2=data.frame(read.delim(paste(protein2_dir, paste(list[i,"uid2"],"_info.csv",sep = ""),sep="/"), sep="\t",header=T, stringsAsFactors = FALSE))
    temp_file2=temp_file2[temp_file2[,"Entry.Name"]==list[i,"uid2"],]
    list[i,"uid2_Pfam"]=temp_file2[1,"InterPro"]
    list[i,"uid2_scl"]=temp_file2[1,5]
    
    system(paste("rm ",  paste(protein2_dir, paste(list[i,"uid2"],"_info.csv",sep = ""), sep='/')))
  
    
    if(strcmpi(list[i,"uid1_Pfam"],list[i,"uid2_Pfam"])==TRUE){ list[i,"shared_domain"]=TRUE}
    
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}






########################################
# calculate proper interface 

#Assign 3d
for(i in 115:nrow(list)){
  #for(i in 1:5){
  tryCatch({
    print(i)
    if(list[i,"diagonal"]=="TRUE"){
    if(list[i,"interface"]!=0){
      
      #create a directory for the every protein
      complex_dir=" "
      complex_dir=paste("/media/hadeer/Elements/human_mito/benckmark/positive/",paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep="_"),sep="")
      dir.create(complex_dir)
      
      complex_pdb=" "
      complex_pdb=paste(complex_dir,"complex_pdb",sep="/")
      dir.create(complex_pdb)
      

      ##################
      
     # pdb_file=data.frame()
     # pdb_file= read.cif( paste(complex_pdb, paste(list[i,"complex_pdb"],"_",list[i,"complex_chain1"],"_",list[i,"complex_chain2"],".pdb",sep = ""), sep = '/'))
      pdb_file=data.frame()
      pdb_file= read.cif(paste(complex_pdb,"/", list[i,"complex_pdb"],".cif",sep=""))
      
      
      if(nrow(pdb_file$atom)==0){next}
      
        chain1=data.frame()
        chain1 = trim.pdb( pdb_file, chain=as.character(list[i,"complex_chain1"]) ,type = 'ATOM')
        #chain1 = trim.pdb( chain1, elety="CA")
        
        #chain1=chain1[as.character(chain1$atom[,"elety"])=="CA",]
        chain2=data.frame()
        chain2 = trim.pdb(pdb_file, chain=as.character(list[i,"complex_chain2"]) ,type = 'ATOM')
        #chain2 = trim.pdb( chain2, elety="CA")
        #chain2=atom.select(chain2, "calpha")
        
        interface=data.frame()
        na_interface=data.frame()
        ca_interface=data.frame()
        cb_interface=data.frame()
        
        
        #calculate the nearest atom interface
        for (j in 1:nrow(chain1$atom)) {
          #for (j in 1:30) {
          Xj= chain1$atom[j,"x"]
          Yj= chain1$atom[j,"y"]
          Zj= chain1$atom[j,"z"]
          
          for (k in 1:nrow(chain2$atom)){  
            Xk= chain2$atom[k,"x"]
            Yk= chain2$atom[k,"y"]
            Zk= chain2$atom[k,"z"]
            
            dista= round( sqrt((Xj-Xk)^2+(Yj-Yk)^2+(Zj-Zk)^2),digits = 4)
            #print(paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"], chain1$atom[j,"elety"] ,chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"], chain2$atom[k,"elety"], dista, sep = " "))
            #allatom_dista=
            
            if(as.character(chain1$atom[j,"elety"])=="CA" & as.character(chain2$atom[k,"elety"])=="CA"& dista <= 15){
              ca_interface_temp=data.frame(test=c(" "))
              #ca_interface_temp[,"V1"]=paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],",",chain1$atom[j,"elety"],chain2$atom[k,"elety"], ",",dista, sep = " ")
              ca_interface_temp[,"test"]=paste( chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"], chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],sep=" ")
              ca_interface_temp[,"ca_ele"]=paste(chain1$atom[j,"elety"],chain2$atom[k,"elety"], sep=" ")
              ca_interface_temp[,"ca_dist"]=paste(dista)
              ca_interface=rbind(ca_interface,ca_interface_temp)
              
              
            }else if(as.character(chain1$atom[j,"elety"])=="CB" & as.character(chain2$atom[k,"elety"])=="CB"& dista <= 15){
              cb_interface_temp=data.frame(test=c(" "))
              #cb_interface_temp[,"V1"]=paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],",",chain1$atom[j,"elety"],chain2$atom[k,"elety"], ",",dista, sep = " ")
              cb_interface_temp[,"test"]=paste( chain1$atom[j,"chain"],chain1$atom[j,"resno"], chain1$atom[j,"resid"], chain2$atom[k,"chain"],chain2$atom[k,"resno"], chain2$atom[k,"resid"],sep=" ")
              cb_interface_temp[,"cb_ele"]=paste(chain1$atom[j,"elety"],chain2$atom[k,"elety"], sep=" ")
              cb_interface_temp[,"cb_dist"]=paste(dista)
              cb_interface=rbind(cb_interface,cb_interface_temp)
              
            }else if(dista <= 5){
              na_interface_temp=data.frame(test=c(" "))
              #na_interface_temp[,"V1"]=paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],",",chain1$atom[j,"elety"],chain2$atom[k,"elety"], ",",dista, sep = " ")
              na_interface_temp[,"test"]=paste( chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],sep=" ")
              na_interface_temp[,"ele"]=paste(chain1$atom[j,"elety"],chain2$atom[k,"elety"], sep=" ")
              na_interface_temp[,"na_dist"]=paste(dista)
              na_interface=rbind(na_interface,na_interface_temp)}
          }
          
        } 
        
        if(nrow(na_interface)>0){ 
          na_interface=na_interface[order(na_interface[,"na_dist"], decreasing=FALSE),]
          na_interface=na_interface[!duplicated(na_interface[,"test"]),]
          interface=merge(na_interface,ca_interface,by.x="test",by.y="test", all.x=TRUE)
          interface=merge(interface,cb_interface,by.x="test",by.y="test", all.x=TRUE)
          interface=separate(interface, test, c("chain1", "resno1", "resid1", "chain2", "resno2", "resid2"), sep=" ")
          interface=separate(interface, ele, c("na_atom1","na_atom2"), sep=" ")
          interface=interface[,-c(10,12)]
          #list[i,"na_contacts"]=nrow(interface)
          write.table(interface,paste(complex_pdb, paste(paste("na_5A_interface",list[i,"complex_pdb"],list[i,"complex_chain1"],list[i,"complex_chain2"],sep='_'),".txt", sep = ''), sep = '/'),sep="\t",row.names=FALSE, quote = FALSE)
   
          
          interface_contact=ggplot()
          interface_contact=interface_contact+ geom_point(data=interface, aes(x=as.numeric(interface[,"resno1"]),y=as.numeric(interface[,"resno2"])),color="darkgrey",size=3)
          interface_contact=interface_contact+ xlim(1,max(as.numeric(interface[,"resno1"]), na.rm = TRUE))+  ylim(1,max(as.numeric(interface[,"resno2"]), na.rm = TRUE))
          interface_contact=interface_contact+  xlab("Residues of protein 1")
          interface_contact=interface_contact+ ylab("Residues of protein 2")
          interface_contact=interface_contact+ theme_bw() + theme( axis.title = element_text(size = 16),  axis.text = element_text(size = 14))
          interface_contact
          ggsave(paste("/media/hadeer/Elements/human_mito/evcomplex/positive/figures/diagonal", paste(list[i,"prefix"], "na_5A_interface.png", sep='_'), sep="/"),  scale = 1, width = NA, height = NA, dpi=300)
          
         }
        
        if(nrow(ca_interface)>0){ 
          ca_interface=ca_interface[order(ca_interface[,"ca_dist"], decreasing=FALSE),]
          ca_interface= ca_interface[!duplicated( ca_interface[,"test"]),]
          ca_interface=separate(ca_interface, test, c("chain1", "resno1", "resid1", "chain2", "resno2", "resid2"), sep=" ")
          ca_interface=separate(ca_interface, ca_ele, c("na_atom1","na_atom2"), sep=" ")
          ca_interface=ca_interface[,-c(10,12)]
          #list[i,"ca_contacts"]=nrow( ca_interface)
          write.table(ca_interface,paste(complex_pdb, paste(paste("ca_15A_interface",list[i,"complex_pdb"],list[i,"complex_chain1"],list[i,"complex_chain2"],sep='_'),".txt", sep = ''), sep = '/'),sep="\t",row.names=FALSE, quote = FALSE)
       
         }
        
      list[i,"na_contact_5A"]=as.numeric(nrow(interface))
      list[i,"ca_contact_15A"]=as.numeric(nrow(ca_interface))
    }
  }
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}

write.csv(list, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d_5840.csv", quote = TRUE)

list2=list[list[,"na_contact_5A"]>=20 ,]
list2=list2[ !is.na(list2[,"prefix"]) ,]
write.csv(list2, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d_307_20contact.csv", quote = TRUE)

list3=list[list[,"na_contact_5A"]>=1 ,]
list3=list3[ !is.na(list3[,"prefix"]) ,]
write.csv(list3, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d_494_1contact.csv", quote = TRUE)


list_backup=list
#remove hopeless cases
list=list[!(list[,"interface"]==0 & list[,"complex_3d"]==1) ,]
#this command decreased the set from 7386 to 6281 
################################################

#comparing strings of Pfam domains to find common domains

for(i in 1:nrow(list)){
  #for(i in 1:10){
  tryCatch({
    print(i)
        #create a directory for the every protein
    protein1_dir=" "
    protein1_dir=paste(working_dir,"dataset/mito_proteome/human",list[i,"uid1"],sep="/")
    #dir.create(protein_dir)
   
    #download protein info from the UniProt
    download.file(paste("https://www.uniprot.org/uniprot/?query=",list[i,"uid1"],"&sort=score&columns=id,database(Pfam)&format=tab",sep=""), destfile = paste(protein1_dir, paste(list[i,"uid1"],"_info.csv",sep = ""), sep='/'))
    
    temp_file1=data.frame()
    temp_file1=data.frame(read.delim(paste(protein1_dir, paste(list[i,"uid1"],"_info.csv",sep = ""),sep="/"), sep=",",header=T, stringsAsFactors = FALSE))
    list[i,"uid1_Pfam"]=temp_file1[1,1]
    system(paste("rm ",  paste(protein1_dir, paste(list[i,"uid1"],"_info.csv",sep = ""), sep='/')))
    
  
    
    #create a directory for the every protein
    protein2_dir=" "
    protein2_dir=paste(working_dir,"dataset/mito_proteome/human",list[i,"uid2"],sep="/")
    #dir.create(protein_dir)
    
    #download protein info from the UniProt
    download.file(paste("https://www.uniprot.org/uniprot/?query=",list[i,"uid2"],"&sort=score&columns=id,database(Pfam)&format=tab",sep=""), destfile = paste(protein2_dir, paste(list[i,"uid2"],"_info.csv",sep = ""), sep='/'))
    
    temp_file2=data.frame()
    temp_file2=data.frame(read.delim(paste(protein2_dir, paste(list[i,"uid2"],"_info.csv",sep = ""),sep="/"), sep=",",header=T, stringsAsFactors = FALSE))
    list[i,"uid2_Pfam"]=temp_file2[1,1]
    system(paste("rm ",  paste(protein2_dir, paste(list[i,"uid2"],"_info.csv",sep = ""), sep='/')))
    
    
    
    if(strcmpi(list[i,"uid1_Pfam"],list[i,"uid2_Pfam"])==TRUE ){ list[i,"shared_domain"]=TRUE}
    
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}

#######################
list0=data.frame(read.delim(paste(working_dir,"about_dataset","mito_proteome/human_mito_structure_ppi_prefix_3d_494_1contact.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
list=list0


for(i in 1:nrow(list)){
  #for(i in 1:10){
  tryCatch({
    print(i)
    
    #create a directory for the every protein
    complex_dir=" "
    complex_dir=paste("/media/hadeer/Elements/human_mito/benckmark/positive/",paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep="_"),sep="")
  
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p00/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p01/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p02/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p03/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p04/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p05/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p06/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p07/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p08/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p09/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p10/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p11/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p12/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p13/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
    system(paste("cp /media/hadeer/Elements/human_mito/evcomplex/sepmito_p14/",paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"), "  ",complex_dir, sep=""))
 
    temp_file=data.frame()
    temp_file=data.frame(read.delim(paste(complex_dir,paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"),sep="/"), header=TRUE,sep=",",as.is=TRUE, colClasses="character", stringsAsFactors = FALSE))
    if(nrow(temp_file)>0){list[i,"ECs_check"]="Complete"}
    
         }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}



for(i in 1:nrow(list)){
  #for(i in 1:10){
  tryCatch({
    print(i)
    
    #create a directory for the every protein
    complex_dir=" "
    complex_dir=paste("/media/hadeer/Elements/human_mito/benckmark/positive/",paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep="_"),sep="")

    temp_file=data.frame()
    #temp_file=data.frame(read.delim(paste(complex_dir,paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"),sep="/"), header=TRUE,sep=",",as.is=TRUE, colClasses="character", stringsAsFactors = FALSE))
    #temp_file=data.frame(read.delim(paste("/media/hadeer/Elements/human_mito/evcomplex/positive",list[i,"prefix"],"couplings", paste(list[i,"prefix"],"CouplingScores_inter.csv",sep="_"),sep="/"), header=TRUE,sep=",",as.is=TRUE, colClasses="character", stringsAsFactors = FALSE))
    temp_file=data.frame(read.delim(paste("/media/hadeer/Elements/human_mito/evcomplex/positive",list[i,"prefix"],"compare", paste(list[i,"prefix"],"CouplingScoresCompared_inter.csv",sep="_"),sep="/"), header=TRUE,sep=",",as.is=TRUE, colClasses="character", stringsAsFactors = FALSE))
    
     if(nrow(temp_file)>0){list[i,"ECs_check"]="Complete"}
    
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}
list[,"ECs_check"]=NA

list3=list[is.na(list[,"ECs_check"]),]
list2=list[is.na(list[,"ECs_check"]),]
list4=list3[!(list3[,"prefix"] %in% list2[,"prefix"]),]

write.csv(list2, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d_5840_sublistredo.csv", quote = TRUE)
/media/hadeer/Elements/human_mito/evcomplex/positive/sepmito_00000007/compare/sepmito_00000007_CouplingScoresCompared_inter.csv

/media/hadeer/Elements/human_mito/evcomplex/positive/sepmito_00000007/couplings/sepmito_00000007_CouplingScores.csv

#1: In file(file, "rt") :
 # cannot open file '/media/hadeer/Elements/human_mito/evcomplex/positive/sepmito_00149028/couplings/sepmito_00149028_CouplingScores_inter.csv': No such file or directory
#2: In file(file, "rt") :
 # cannot open file '/media/hadeer/Elements/human_mito/evcomplex/positive/sepmito_00001665/couplings/sepmito_00001665_CouplingScores_inter.csv': No such file or directory

####################################
###ecs of postives 
#####################################

list2=list[list[,"na_contact_5A"]>=20 ,]
list2=list2[ !is.na(list2[,"prefix"]) ,]


for (i in 1:nrow(list2)){
  tryCatch({
    print(i)

    complex_dir=" "
    complex_dir=paste("/media/hadeer/Elements/human_mito/benckmark/positive/",paste(list2[i,"prefix"],sub("\\_.*", "",list2[i,"uid1"]),sub("\\_.*", "",list2[i,"uid2"]),sep="_"),sep="")
    
    ec_list=data.frame()
    ec_list=data.frame(read.delim(paste(complex_dir,paste(list2[i,"prefix"],"CouplingScores_inter.csv",sep="_"),sep="/"), header=TRUE,sep=",",as.is=TRUE, colClasses="character", stringsAsFactors = FALSE))
    
    ec_list= ec_list[order( as.numeric(ec_list[,"cn"]), decreasing = TRUE),]
    
    list2[i,"kurtosis"]=kurtosis(as.numeric(ec_list[,"cn"]))
    list2[i,"skewness"]=skewness(as.numeric(ec_list[,"cn"]))
    list2[i,"squared_skewness"]=( list2[i,"skewness"])^2
    list2[i,"max_ecs"]=max(as.numeric(ec_list[,"cn"]))
    
    
  },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}


plot= ggplot() + ggtitle("")
plot= plot+ geom_density(list2, aes(x=as.numeric( list2[,"skewness"])),color="lightblue")
plot 


plot= ggplot() + ggtitle("")
plot= plot+ geom_point(data=list2,aes(x=as.numeric( list2[,"skewness"]), y=as.numeric(list2[,"kurtosis"])), show.legend = TRUE)
#plot= plot+ geom_point(data=list,aes(x=as.numeric( list[,"skewness"]), y=as.numeric(list[,"kurtosis"]),color="TP_PPI"),color="darkslateblue", show.legend = TRUE)
plot= plot+ geom_hline(yintercept=10, linetype="dashed", color = "red", size=0.5)
plot= plot+ xlim(-2.5,5)+ylim(-2.5,20)
plot= plot+ geom_point(aes(x=0, y=0), colour="black",size=4,)
plot= plot+ xlab("Skewness") +ylab("Excess Kurtosis") +theme_bw()
plot= plot+ theme(axis.text=element_text(size=10),axis.title=element_text(size=16))
ggMarginal(plot,type='histogram', color="black",fill="deepskyblue4",alpha=0.7)
plot

#ggsave("EkurtosisVSskewness_p3_density.png", plot = last_plot(), path = "/Hadeer/hybrid_Xlmsecs/stat_plot", scale = 1, width = NA, height = NA, dpi=300)

#I ran until here
##############################################################
















#####
# I will try to check the second option in blast results 

for(i in 1:nrow(list)){
  #for(i in 1:10){
  tryCatch({
    print(i)
    
    if( is.na(list[i,"interface"]) || list[i,"interface"]==0){
      #if( is.na(list[i,"interface"])){
      
      #check the second unique pdb in the blast search  
      temp_file=data.frame()
      temp_file=data.frame(read.delim(paste("/Hadeer/hybrid_Xlmsecs/dataset/mito_complex/struc_PPIs/human/",list[i,"interactorA_uid"],"-",list[i,"interactorB_uid"],".csv"), header=TRUE,sep=",", stringsAsFactors = FALSE))
      temp_file=temp_file[!duplicated(temp_file[,"pdb"]),]
      temp_file1=temp_file[temp_file[,"pdb"]==temp_file[2,"pdb"],]
      
      temp_file=data.frame()
      temp_file=data.frame(read.delim(paste("/Hadeer/hybrid_Xlmsecs/dataset/mito_complex/struc_PPIs/human/",list[i,"interactorA_uid"],"-",list[i,"interactorB_uid"],".csv"), header=TRUE,sep=",", stringsAsFactors = FALSE))
      temp_file=temp_file[temp_file[,"pdb"]==temp_file1[1,"pdb"],]
      
      download.file(paste(paste("https://files.rcsb.org/view", temp_file[1,"pdb"], sep = '/'), "pdb", sep = '.'), destfile = paste("/Hadeer/hybrid_Xlmsecs/dataset/mito_complex/struc_PPIs/human/", temp_file[1,"pdb"],".pdb",sep=""))
      pdb_file=data.frame()
      pdb_file= read.pdb(paste("/Hadeer/hybrid_Xlmsecs/dataset/mito_complex/struc_PPIs/human/", temp_file[1,"pdb"],".pdb",sep=""))
      if(nrow(pdb_file$atom)==0){next}
      
      for(l in 1:nrow(temp_file)){
        chain1=data.frame()
        chain1 = trim.pdb( pdb_file, chain=as.character(temp_file[l,"chain.x"]) ,type = 'ATOM')
        chain1 = trim.pdb( chain1, elety="CA")
        
        #chain1=chain1[as.character(chain1$atom[,"elety"])=="CA",]
        chain2=data.frame()
        chain2 = trim.pdb(pdb_file, chain=as.character(temp_file[l,"chain.y"]) ,type = 'ATOM')
        chain2 = trim.pdb( chain2, elety="CA")
        #chain2=atom.select(chain2, "calpha")
        
        na_interface=data.frame()
        #calculate the nearest atom interface
        for (j in 1:nrow(chain1$atom)) {
          #for (j in 1:30) {
          Xj= chain1$atom[j,"x"]
          Yj= chain1$atom[j,"y"]
          Zj= chain1$atom[j,"z"]
          
          for (k in 1:nrow(chain2$atom)){  
            Xk= chain2$atom[k,"x"]
            Yk= chain2$atom[k,"y"]
            Zk= chain2$atom[k,"z"]
            
            dista= round( sqrt((Xj-Xk)^2+(Yj-Yk)^2+(Zj-Zk)^2),digits = 4)
            
            if(dista <= 10){
              na_interface_temp=data.frame(test=c(" "))
              na_interface_temp[,"test"]=paste( chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],sep=" ")
              na_interface_temp[,"ele"]=paste(chain1$atom[j,"elety"],chain2$atom[k,"elety"], sep=" ")
              na_interface_temp[,"na_dist"]=paste(dista)
              na_interface=rbind(na_interface,na_interface_temp)}
          }
          
          
          if(nrow(na_interface)>10){break}
        }
        
        if(nrow(na_interface)>0){ temp_file[l,"interface"]=TRUE}else{temp_file[l,"interface"]=FALSE}
      }
      
      
      list[i,"interface"]=as.numeric(nrow(temp_file[temp_file[,"interface"]==TRUE,]))
      
      temp_file1=data.frame          
      temp_file1=temp_file[temp_file[,"interface"]==TRUE & !duplicated(temp_file[,"chain.x"]),]
      temp_file2=data.frame 
      temp_file2=temp_file[temp_file[,"interface"]==TRUE & !duplicated(temp_file[,"chain.y"]),]  
      temp_file[1,"interface_architect"]=paste(nrow(temp_file1), ":" ,nrow(temp_file2), sep=" ")
      
      
      list[i,"complex_psb"]=temp_file[1,"pdb"]
      list[i,"complex_chain1"]=temp_file[1,"chain.x"]
      list[i,"complex_chain2"]=temp_file[1,"chain.y"]
      list[i,"no_pairs"]=nrow(temp_file)
      list[i,"interface_architect"]=temp_file[1,"interface_architect"]
      
      system(paste("rm ", paste("/Hadeer/hybrid_Xlmsecs/dataset/mito_complex/struc_PPIs/human/", temp_file[1,"pdb"],".pdb",sep="") ))
      
    }
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}






#remove hopeless cases
list=list[!(list[,"interface"]==0 & list[,"complex_3d"]==2) ,]
list=list[!is.na(list[,"interactorA_uid"]) ,]
#this command decreased the set from  6281 to 4948

list=list[!(list[,"interface"]==0) ,]
list=list[!is.na(list[,"interactorA_uid"]) ,]
#this command decreased the set from   4948 to 3688
write.csv(list_backup, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_positive_mito_ppi_raw.csv")

write.csv(list, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_positive_mito_ppi.csv")



list=list[,-c(1:17,20,21,25,35,47)]
list_2=list[list[,"na_contact_5A"]>=1,]
names(list)[7]="uid1"
names(list)[15]="uid2"


##############################################################
write.csv(list, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d.csv", quote = TRUE)

list0=data.frame(read.delim(paste(working_dir,"about_dataset","mito_proteome/human_mito_structure_ppi_prefix_3d.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
list=list0

#Assign 3d
for(i in 1:nrow(list)){
  #for(i in 1:5){
  tryCatch({
    print(i)
    
    if( is.na(list[i,"interface"]) || list[i,"interface"]==0){
      
    
    #create a directory for the every protein
    protein1_dir=" "
    protein1_dir=paste(working_dir,"dataset/mito_proteome/human",list[i,"uid1"],sep="/")
    #dir.create(protein_dir)
    
    temp_file1=data.frame()
    temp_file1=data.frame(read.delim(paste(protein1_dir,paste(list[i,"uid1"], "blastp.out",sep="_"), sep='/'), header=FALSE, stringsAsFactors = FALSE))
    names(temp_file1)[1:12]= c("query", "subject", "%identity", "alignment length", "mismatches", "gap opens", "q. start:", "q. end", "s. start", "s. end", "evalue", "bit score")
    temp_file1[,"pdb"]=gsub("_.*", "\\4", temp_file1[,"subject"])
    temp_file1[,"chain"]=gsub(".*_", "\\1", temp_file1[,"subject"])
    temp_file1=temp_file1[!duplicated(temp_file1[,"subject"]),]
    temp_file1=temp_file1[temp_file1[,"%identity"]>=50,]
    
    #create a directory for the every protein
    protein2_dir=" "
    protein2_dir=paste(working_dir,"dataset/mito_proteome/human",list[i,"uid2"],sep="/")
    #dir.create(protein_dir)
    
    temp_file2=data.frame()
    temp_file2=data.frame(read.delim(paste(protein2_dir,paste(list[i,"uid2"], "blastp.out",sep="_"), sep='/'), header=FALSE, stringsAsFactors = FALSE))
    names(temp_file2)[1:12]= c("query", "subject", "%identity", "alignment length", "mismatches", "gap opens", "q. start:", "q. end", "s. start", "s. end", "evalue", "bit score")
    temp_file2[,"pdb"]=gsub("_.*", "\\4", temp_file2[,"subject"])
    temp_file2[,"chain"]=gsub(".*_", "\\1", temp_file2[,"subject"])
    temp_file2=temp_file2[!duplicated(temp_file2[,"subject"]),]
    temp_file2=temp_file2[temp_file2[,"%identity"]>=50,]
    
    #remove case where two proteins wee mapped to the same chain
    temp_file=data.frame()
    temp_file=merge(temp_file1,temp_file2,by.x="pdb", by.y="pdb")
    temp_file= temp_file[temp_file[,"subject.x"]!=temp_file[,"subject.y"],]   
    
    #remove cases where  two proteins wee mapped to the same chain intercahnagably
    temp_file[,"test1"]=paste(temp_file[,"subject.x"], temp_file[,"subject.y"],sep=" ")
    temp_file[,"test2"]=paste(temp_file[,"subject.y"], temp_file[,"subject.x"],sep=" ")
    
    temp_file[,"delete"]=" "
    for(k in 1:nrow(temp_file)){
      for(l in 1:nrow(temp_file)){
        if(temp_file[k,"test1"]==temp_file[l,"test2"]){
          temp_file[k,"delete"]=TRUE
          temp_file[l,"delete"]=TRUE
          print("case detected")
        }
      }
    } 
    
    temp_file= temp_file[temp_file[,"delete"]!=TRUE,]   
    # temp_file= temp_file[ order(as.numeric(temp_file[,"evalue.x"]),as.numeric(temp_file[,"evalue.y"])),]    
    temp_file= temp_file[ order(as.numeric(temp_file[,"bit score.x"]),as.numeric(temp_file[,"bit score.y"]),decreasing =TRUE ),]    
    
    #create a directory for the every protein
    complex_dir=" "
    complex_dir=paste("/media/hadeer/Elements/human_mito/benckmark/positive/",paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep="_"),sep="")
    dir.create(complex_dir)
    
    complex_pdb=" "
    complex_pdb=paste(complex_dir,"complex_pdb",sep="/")
    dir.create(complex_pdb)
    
    if(nrow(temp_file)!=0){
      
      write.csv( temp_file, file =paste(complex_pdb,paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),"blastout.csv",sep="_"),sep="/"))
    }
    
    ##################
    
    temp_file=data.frame()
    temp_file=data.frame(read.delim(paste(complex_pdb,paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),"blastout.csv",sep="_"),sep="/"), header=TRUE,sep=",",as.is=TRUE, colClasses="character", stringsAsFactors = FALSE))
    temp_file=temp_file[temp_file[,"pdb"]==temp_file[1,"pdb"],]
    download.file(paste(paste("https://files.rcsb.org/view", temp_file[1,"pdb"], sep = '/'), "pdb", sep = '.'), destfile = paste( complex_pdb,"/", temp_file[1,"pdb"],".pdb",sep=""))
    pdb_file=data.frame()
    pdb_file= read.pdb(paste(complex_pdb,"/", temp_file[1,"pdb"],".pdb",sep=""))
    if(nrow(pdb_file$atom)==0){next}
    
    for(l in 1:nrow(temp_file)){
      chain1=data.frame()
      chain1 = trim.pdb( pdb_file, chain=as.character(temp_file[l,"chain.x"]) ,type = 'ATOM')
      #chain1 = trim.pdb( chain1, elety="CA")
      
      #chain1=chain1[as.character(chain1$atom[,"elety"])=="CA",]
      chain2=data.frame()
      chain2 = trim.pdb(pdb_file, chain=as.character(temp_file[l,"chain.y"]) ,type = 'ATOM')
      #chain2 = trim.pdb( chain2, elety="CA")
      #chain2=atom.select(chain2, "calpha")
      
      interface=data.frame()
      na_interface=data.frame()
      ca_interface=data.frame()
      cb_interface=data.frame()
      
      
      #calculate the nearest atom interface
      for (j in 1:nrow(chain1$atom)) {
        #for (j in 1:30) {
        Xj= chain1$atom[j,"x"]
        Yj= chain1$atom[j,"y"]
        Zj= chain1$atom[j,"z"]
        
        for (k in 1:nrow(chain2$atom)){  
          Xk= chain2$atom[k,"x"]
          Yk= chain2$atom[k,"y"]
          Zk= chain2$atom[k,"z"]
          
          dista= round( sqrt((Xj-Xk)^2+(Yj-Yk)^2+(Zj-Zk)^2),digits = 4)
          #print(paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"], chain1$atom[j,"elety"] ,chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"], chain2$atom[k,"elety"], dista, sep = " "))
          #allatom_dista=
          
          if(as.character(chain1$atom[j,"elety"])=="CA" & as.character(chain2$atom[k,"elety"])=="CA"& dista <= 15){
            ca_interface_temp=data.frame(test=c(" "))
            #ca_interface_temp[,"V1"]=paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],",",chain1$atom[j,"elety"],chain2$atom[k,"elety"], ",",dista, sep = " ")
            ca_interface_temp[,"test"]=paste( chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"], chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],sep=" ")
            ca_interface_temp[,"ca_ele"]=paste(chain1$atom[j,"elety"],chain2$atom[k,"elety"], sep=" ")
            ca_interface_temp[,"ca_dist"]=paste(dista)
            ca_interface=rbind(ca_interface,ca_interface_temp)
            
            
          }else if(as.character(chain1$atom[j,"elety"])=="CB" & as.character(chain2$atom[k,"elety"])=="CB"& dista <= 15){
            cb_interface_temp=data.frame(test=c(" "))
            #cb_interface_temp[,"V1"]=paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],",",chain1$atom[j,"elety"],chain2$atom[k,"elety"], ",",dista, sep = " ")
            cb_interface_temp[,"test"]=paste( chain1$atom[j,"chain"],chain1$atom[j,"resno"], chain1$atom[j,"resid"], chain2$atom[k,"chain"],chain2$atom[k,"resno"], chain2$atom[k,"resid"],sep=" ")
            cb_interface_temp[,"cb_ele"]=paste(chain1$atom[j,"elety"],chain2$atom[k,"elety"], sep=" ")
            cb_interface_temp[,"cb_dist"]=paste(dista)
            cb_interface=rbind(cb_interface,cb_interface_temp)
            
          }else if(dista <= 5){
            na_interface_temp=data.frame(test=c(" "))
            #na_interface_temp[,"V1"]=paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],",",chain1$atom[j,"elety"],chain2$atom[k,"elety"], ",",dista, sep = " ")
            na_interface_temp[,"test"]=paste( chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],sep=" ")
            na_interface_temp[,"ele"]=paste(chain1$atom[j,"elety"],chain2$atom[k,"elety"], sep=" ")
            na_interface_temp[,"na_dist"]=paste(dista)
            na_interface=rbind(na_interface,na_interface_temp)}
        }
        
      } 
      
      if(nrow(na_interface)>0){ 
        na_interface=na_interface[order(na_interface[,"na_dist"], decreasing=FALSE),]
        na_interface=na_interface[!duplicated(na_interface[,"test"]),]
        interface=merge(na_interface,ca_interface,by.x="test",by.y="test", all.x=TRUE)
        interface=merge(interface,cb_interface,by.x="test",by.y="test", all.x=TRUE)
        interface=separate(interface, test, c("chain1", "resno1", "resid1", "chain2", "resno2", "resid2"), sep=" ")
        interface=separate(interface, ele, c("na_atom1","na_atom2"), sep=" ")
        interface=interface[,-c(10,12)]
        #list[i,"na_contacts"]=nrow(interface)
        write.table(interface,paste(complex_pdb, paste(paste("na_5A_interface",temp_file[l,"pdb"],temp_file[l,"chain.x"],temp_file[l,"chain.y"],sep='_'),".txt", sep = ''), sep = '/'),sep="\t",row.names=FALSE, quote = FALSE)
      }
      
      
      if(nrow(ca_interface)>0){ 
        ca_interface=ca_interface[order(ca_interface[,"ca_dist"], decreasing=FALSE),]
        ca_interface= ca_interface[!duplicated( ca_interface[,"test"]),]
        ca_interface=separate(ca_interface, test, c("chain1", "resno1", "resid1", "chain2", "resno2", "resid2"), sep=" ")
        ca_interface=separate(ca_interface, ca_ele, c("na_atom1","na_atom2"), sep=" ")
        ca_interface=ca_interface[,-c(10,12)]
        #list[i,"ca_contacts"]=nrow( ca_interface)
        write.table(ca_interface,paste(complex_pdb, paste(paste("ca_15A_interface",temp_file[l,"pdb"],temp_file[l,"chain.x"],temp_file[l,"chain.y"],sep='_'),".txt", sep = ''), sep = '/'),sep="\t",row.names=FALSE, quote = FALSE)
      }
      
      
      temp_file[l,"na_contacts"]=as.numeric(nrow(interface))
      temp_file[l,"ca_contacts"]=as.numeric(nrow(ca_interface))
      
      if(nrow(na_interface)>0){ temp_file[l,"interface"]=TRUE}else{temp_file[l,"interface"]=FALSE}
      
    }
    
    
    
    list[i,"interface"]=as.numeric(nrow(temp_file[temp_file[,"interface"]==TRUE,]))
    
    temp_file1=data.frame          
    temp_file1=temp_file[temp_file[,"interface"]==TRUE & !duplicated(temp_file[,"chain.x"]),]
    temp_file2=data.frame 
    temp_file2=temp_file[temp_file[,"interface"]==TRUE & !duplicated(temp_file[,"chain.y"]),]  
    temp_file[1,"interface_architect"]=paste(nrow(temp_file1), ":" ,nrow(temp_file2), sep=" ")
    
    
    list[i,"complex_psb"]=temp_file[1,"pdb"]
    list[i,"complex_chain1"]=temp_file[1,"chain.x"]
    list[i,"complex_chain2"]=temp_file[1,"chain.y"]
    list[i,"no_pairs"]=nrow(temp_file)
    list[i,"interface_architect"]=temp_file[1,"interface_architect"]
    list[i,"na_contact_5A"]=temp_file[1,"na_contacts"]
    list[i,"ca_contact_15A"]=temp_file[1,"ca_contacts"]
    }
    
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}


for(i in 1:nrow(list)){
  #for(i in 1:5){
  tryCatch({
    print(i)
    
    complex_dir=" "
    complex_dir=paste(working_dir,"/dataset/mito_complex/human/benchmark/",paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep="_"),sep="")
    dir.create(complex_dir)
    
    complex_pdb=" "
    complex_pdb=paste(complex_dir,"complex_pdb",sep="/")
    
    temp_file=data.frame()
    temp_file=data.frame(read.delim(paste(complex_pdb,paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),"blastout.csv",sep="_"),sep="/"), header=TRUE,sep=",",as.is=TRUE, colClasses="character", stringsAsFactors = FALSE))
    
    
    list[i,"n_complex"]=as.numeric(nrow(temp_file))
    
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}



list=list[!is.na(list[,"prefix"]),]
list=list[order(list[,"prefix"]),]

write.csv(list, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d.csv", quote = TRUE)



for(i in 1:nrow(list)){
  #for(i in 1:5){
  tryCatch({
    print(i)
    complex_dir=" "
    complex_dir=paste(working_dir,"/dataset/mito_complex/human/benchmark/",paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep="_"),sep="")
    dir.create(complex_dir)
    
    complex_pdb=" "
    complex_pdb=paste(complex_dir,"complex_pdb",sep="/")
    
    evcomplex=" "
    evcomplex=paste(complex_dir,"evcomplex",sep="/")
    dir.create(evcomplex)
    
    system(paste("cp ", paste("/media/hadeer/Elements/human_mito/evcomplex/sepmito_p13/", list[i,"prefix"],"_CouplingScores_inter.csv",sep=""), evcomplex,sep=" " ))
    
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}


##############################################################

#Assign 3d for the second entry 
for(i in 1:nrow(list)){
  #  for(i in 1:10){
  tryCatch({
    print(i)
    
    
    if(list[i,"n_complex"]>=2 & list[i,"interface"]==0 ){
      
      #create a directory for the every protein
      complex_dir=" "
      complex_dir=paste(working_dir,"/dataset/mito_complex/human/benchmark/",paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep="_"),sep="")
      dir.create(complex_dir)
      
      complex_pdb=" "
      complex_pdb=paste(complex_dir,"complex_pdb",sep="/")
      
      ##################
      
      temp_file=data.frame()
      temp_file=data.frame(read.delim(paste(complex_pdb,paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),"blastout.csv",sep="_"),sep="/"), header=TRUE,sep=",",as.is=TRUE, colClasses="character", stringsAsFactors = FALSE))
      temp_file=temp_file[temp_file[,"pdb"]!=temp_file[1,"pdb"],]
      if(nrow(temp_file)==0){next}
      download.file(paste(paste("https://files.rcsb.org/view", temp_file[1,"pdb"], sep = '/'), "pdb", sep = '.'), destfile = paste( complex_pdb,"/", temp_file[1,"pdb"],".pdb",sep=""))
      pdb_file=data.frame()
      pdb_file= read.pdb(paste(complex_pdb,"/", temp_file[1,"pdb"],".pdb",sep=""))
      if(nrow(pdb_file$atom)==0){next}
      
      for(l in 1:nrow(temp_file)){
        chain1=data.frame()
        chain1 = trim.pdb( pdb_file, chain=as.character(temp_file[l,"chain.x"]) ,type = 'ATOM')
        #chain1 = trim.pdb( chain1, elety="CA")
        
        #chain1=chain1[as.character(chain1$atom[,"elety"])=="CA",]
        chain2=data.frame()
        chain2 = trim.pdb(pdb_file, chain=as.character(temp_file[l,"chain.y"]) ,type = 'ATOM')
        #chain2 = trim.pdb( chain2, elety="CA")
        #chain2=atom.select(chain2, "calpha")
        
        interface=data.frame()
        na_interface=data.frame()
        ca_interface=data.frame()
        cb_interface=data.frame()
        
        
        #calculate the nearest atom interface
        for (j in 1:nrow(chain1$atom)) {
          #for (j in 1:30) {
          Xj= chain1$atom[j,"x"]
          Yj= chain1$atom[j,"y"]
          Zj= chain1$atom[j,"z"]
          
          for (k in 1:nrow(chain2$atom)){  
            Xk= chain2$atom[k,"x"]
            Yk= chain2$atom[k,"y"]
            Zk= chain2$atom[k,"z"]
            
            dista= round( sqrt((Xj-Xk)^2+(Yj-Yk)^2+(Zj-Zk)^2),digits = 4)
            #print(paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"], chain1$atom[j,"elety"] ,chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"], chain2$atom[k,"elety"], dista, sep = " "))
            #allatom_dista=
            
            if(as.character(chain1$atom[j,"elety"])=="CA" & as.character(chain2$atom[k,"elety"])=="CA"& dista <= 15){
              ca_interface_temp=data.frame(test=c(" "))
              #ca_interface_temp[,"V1"]=paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],",",chain1$atom[j,"elety"],chain2$atom[k,"elety"], ",",dista, sep = " ")
              ca_interface_temp[,"test"]=paste( chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"], chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],sep=" ")
              ca_interface_temp[,"ca_ele"]=paste(chain1$atom[j,"elety"],chain2$atom[k,"elety"], sep=" ")
              ca_interface_temp[,"ca_dist"]=paste(dista)
              ca_interface=rbind(ca_interface,ca_interface_temp)
              
              
            }else if(as.character(chain1$atom[j,"elety"])=="CB" & as.character(chain2$atom[k,"elety"])=="CB"& dista <= 15){
              cb_interface_temp=data.frame(test=c(" "))
              #cb_interface_temp[,"V1"]=paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],",",chain1$atom[j,"elety"],chain2$atom[k,"elety"], ",",dista, sep = " ")
              cb_interface_temp[,"test"]=paste( chain1$atom[j,"chain"],chain1$atom[j,"resno"], chain1$atom[j,"resid"], chain2$atom[k,"chain"],chain2$atom[k,"resno"], chain2$atom[k,"resid"],sep=" ")
              cb_interface_temp[,"cb_ele"]=paste(chain1$atom[j,"elety"],chain2$atom[k,"elety"], sep=" ")
              cb_interface_temp[,"cb_dist"]=paste(dista)
              cb_interface=rbind(cb_interface,cb_interface_temp)
              
            }else if(dista <= 5){
              na_interface_temp=data.frame(test=c(" "))
              #na_interface_temp[,"V1"]=paste(chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],",",chain1$atom[j,"elety"],chain2$atom[k,"elety"], ",",dista, sep = " ")
              na_interface_temp[,"test"]=paste( chain1$atom[j,"chain"], chain1$atom[j,"resno"], chain1$atom[j,"resid"],chain2$atom[k,"chain"], chain2$atom[k,"resno"], chain2$atom[k,"resid"],sep=" ")
              na_interface_temp[,"ele"]=paste(chain1$atom[j,"elety"],chain2$atom[k,"elety"], sep=" ")
              na_interface_temp[,"na_dist"]=paste(dista)
              na_interface=rbind(na_interface,na_interface_temp)}
          }
          
        } 
        
        if(nrow(na_interface)>0){ 
          na_interface=na_interface[order(na_interface[,"na_dist"], decreasing=FALSE),]
          na_interface=na_interface[!duplicated(na_interface[,"test"]),]
          interface=merge(na_interface,ca_interface,by.x="test",by.y="test", all.x=TRUE)
          interface=merge(interface,cb_interface,by.x="test",by.y="test", all.x=TRUE)
          interface=separate(interface, test, c("chain1", "resno1", "resid1", "chain2", "resno2", "resid2"), sep=" ")
          interface=separate(interface, ele, c("na_atom1","na_atom2"), sep=" ")
          interface=interface[,-c(10,12)]
          #list[i,"na_contacts"]=nrow(interface)
          write.table(interface,paste(complex_pdb, paste(paste("na_5A_interface",temp_file[l,"pdb"],temp_file[l,"chain.x"],temp_file[l,"chain.y"],sep='_'),".txt", sep = ''), sep = '/'),sep="\t",row.names=FALSE, quote = FALSE)
        }
        
        
        if(nrow(ca_interface)>0){ 
          ca_interface=ca_interface[order(ca_interface[,"ca_dist"], decreasing=FALSE),]
          ca_interface= ca_interface[!duplicated( ca_interface[,"test"]),]
          ca_interface=separate(ca_interface, test, c("chain1", "resno1", "resid1", "chain2", "resno2", "resid2"), sep=" ")
          ca_interface=separate(ca_interface, ca_ele, c("na_atom1","na_atom2"), sep=" ")
          ca_interface=ca_interface[,-c(10,12)]
          #list[i,"ca_contacts"]=nrow( ca_interface)
          write.table(ca_interface,paste(complex_pdb, paste(paste("ca_15A_interface",temp_file[l,"pdb"],temp_file[l,"chain.x"],temp_file[l,"chain.y"],sep='_'),".txt", sep = ''), sep = '/'),sep="\t",row.names=FALSE, quote = FALSE)
        }
        
        
        temp_file[l,"na_contacts"]=as.numeric(nrow(interface))
        temp_file[l,"ca_contacts"]=as.numeric(nrow(ca_interface))
        
        if(nrow(na_interface)>0){ temp_file[l,"interface"]=TRUE}else{
          temp_file[l,"interface"]=FALSE
          system(paste("rm ", paste( complex_pdb,"/", temp_file[1,"pdb"],".pdb",sep="") ,sep=" " ))
        }
        
      }
      
    }
    
    if(nrow(na_interface)>0){ 
      print(paste(i,list[i,"prefix"],sep=" "))
      
      list[i,"interface"]=as.numeric(nrow(temp_file[temp_file[,"interface"]==TRUE,]))
      
      temp_file1=data.frame          
      temp_file1=temp_file[temp_file[,"interface"]==TRUE & !duplicated(temp_file[,"chain.x"]),]
      temp_file2=data.frame 
      temp_file2=temp_file[temp_file[,"interface"]==TRUE & !duplicated(temp_file[,"chain.y"]),]  
      temp_file[1,"interface_architect"]=paste(nrow(temp_file1), ":" ,nrow(temp_file2), sep=" ")
      
      
      list[i,"complex_psb"]=temp_file[1,"pdb"]
      list[i,"complex_chain1"]=temp_file[1,"chain.x"]
      list[i,"complex_chain2"]=temp_file[1,"chain.y"]
      list[i,"no_pairs"]=nrow(temp_file)
      list[i,"interface_architect"]=temp_file[1,"interface_architect"]
      list[i,"na_contact_5A"]=temp_file[1,"na_contacts"]
      list[i,"ca_contact_15A"]=temp_file[1,"ca_contacts"]
      
    }else{
      next}
    
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}





###################################################################




##############################################################

#comparing strings of Pfam domains to find common domains
struc_ppi=list
struc_ppi[,"shared_domain"]=" "
for(i in 1:nrow(struc_ppi)){
  #for(i in 1:10){
  tryCatch({
    print(i)
    if(struc_ppi[i,"interactorA_Pfam_domains"]!=""|| struc_ppi[i,"interactorB_Pfam_domains"]!=""){
      if(strcmp(struc_ppi[i,"interactorA_Pfam_domains"],struc_ppi[i,"interactorB_Pfam_domains"])==TRUE ){ struc_ppi[i,"shared_domain"]=TRUE}
    }
    
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
}

struc_ppi=struc_ppi[struc_ppi[,"shared_domain"]!=TRUE,]
#this command decreased the set from  3688 to 3302

write.csv(struc_ppi, file ="/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_positive_mito_ppi_selected.csv")

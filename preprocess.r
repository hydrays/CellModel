require(XML)
Cfile <- file("commands", "w")

cmd1 <- paste("cp config.xml config.old")
system(cmd1)

##output_file_prefix <- format(Sys.time(), "%Y%m%d%H%M")
output_file_prefix <- "braf8"
output_path = "out"
rep = 10

output_dir = "out"
cell_position_file = "cellPos.txt"
Nstep = 2000
dt = 0.01
eta = 20
q0 = 1.53
a = 1.2
b = 0.8
a_braf = 0.5
b_braf = 0.5
braf_mosaic_percentage = 0.4
mu = 0.05
stretching_force = 2.0
H_rate = 0.004
flag_output_cell = 1
flag_output_force = 0
flag_record_final_state = 0

##for ( stretching_force in seq(5.5, 15, by=0.5) )
for ( i in seq(10) )
{
    rep = rep + 1
    cat("\n ***********run********* ", rep, "\n")
    FolderName <- paste("case_", output_file_prefix, "_", as.character(rep), sep='')
    dir.create(FolderName)

    ## The following part is to produce "config.xml"
    d <- xmlTreeParse("config.xml")
    r <- xmlRoot(d)
    xmlValue(r[["output_dir"]]) <- output_dir
    xmlValue(r[["cell_position_file"]]) <- cell_position_file
    xmlValue(r[["Nstep"]]) <- Nstep
    xmlValue(r[["dt"]]) <- dt
    xmlValue(r[["eta"]]) <- eta
    xmlValue(r[["q0"]]) <- q0
    xmlValue(r[["a"]]) <- a
    xmlValue(r[["b"]]) <- b
    xmlValue(r[["a_braf"]]) <- a_braf
    xmlValue(r[["b_braf"]]) <- b_braf
    xmlValue(r[["braf_mosaic_percentage"]]) <- braf_mosaic_percentage
    xmlValue(r[["mu"]]) <- mu
    xmlValue(r[["stretching_force"]]) <- stretching_force
    xmlValue(r[["H_rate"]]) <- H_rate
    xmlValue(r[["flag_output_cell"]]) <- flag_output_cell
    xmlValue(r[["flag_output_force"]]) <- flag_output_force
    xmlValue(r[["flag_record_final_state"]]) <- flag_record_final_state
    saveXML(r, "config.xml")

    ## Copy the seqfile.txt into folder
    cmd2 <- paste("cp config.xml", FolderName)
    cmd3 <- paste("mkdir -p ", FolderName, "/out", sep='')
    cmd4 <- paste("cp PhysModel", FolderName)
    cmd5 <- paste("cp cellPos.txt", FolderName)
    system(cmd2)
    system(cmd3)
    system(cmd4)
    system(cmd5)

    ##rtime <- runif(1, 20, 200)
    ##runcmd <- paste("sleep ", rtime, "; cd", FolderName, "; screen PhysModel; cd ..")
    runcmd <- paste("cd", FolderName, "; screen -d -m PhysModel; cd ..")
    ##cat(paste("sleep ", rtime, "; cd", FolderName, "; PhysModel; cd .."), file="run")
    ##runcmd <- paste("screen -d -m ", FolderName, "/PhysModel", sep='')
    ##runcmd <- paste("bsub ", FolderName, "/run -o output_", i, sep='')	
    ##runcmd <- paste("cd", FolderName, "; pwd; BigCAT; cd ..;")
    writeLines(runcmd, Cfile)

}

cmd1 <- paste("mv config.old config.xml")
system(cmd1)

close(Cfile)


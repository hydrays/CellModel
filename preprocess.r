require(XML)
Cfile <- file("commands", "w")

cmd1 <- paste("cp config.xml config.old")
system(cmd1)

##output_file_prefix <- format(Sys.time(), "%Y%m%d%H%M")
output_file_prefix <- "basemodel2_a135_showcell2"
output_path = "out"
rep = 0

output_dir = "out"
cell_position_file = "cellPos.txt"
Nstep = 1600
dt = 0.01
eta = 20
two_population_model = 1
q0= 1.5
a = 1.35
b = 1.0
a_braf = 0.5
b_braf = 0.5
braf_mosaic_percentage = 0.0
mu = 0.01
stretching_force = 3.0
H_rate = 0.004
flag_output_cell = 1
flag_output_force = 0
flag_record_final_state = 0

stretching_force_list <- c(1, 3, 7, 10, 14, 16)
##for ( stretching_force in stretching_force_list )
##for ( stretching_force in c(seq(15, 15, length.out=10), seq(10, 10, length.out=10) ))
##for ( i in seq(10) )
{
    ##stretching_force = 4
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
    xmlValue(r[["two_population_model"]]) <- two_population_model
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

    ## screen, on PC
    runcmd <- paste("cd", FolderName, "; screen -d -m PhysModel; cd ..")
    ## runcmd <- paste("screen -d -m ", FolderName, "/PhysModel", sep='')

    ## ## bsub, on server
    ## runcmd <- paste("bsub ", FolderName, "/PhysModel -o output_", rep, sep='')	
    ## runcmd <- paste("cd ", FolderName, "; bsub PhysModel -o output_", rep, "; cd ..", sep='')	

    writeLines(runcmd, Cfile)

}

cmd1 <- paste("mv config.old config.xml")
system(cmd1)

close(Cfile)


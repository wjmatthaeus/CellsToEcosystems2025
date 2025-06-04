##data inputs generated using point simulation with Paleo-BGC

pkgs <- c("readr", "ggplot2", "dplyr", "readxl",
          "tibble", "tidyr", "stringr", "raster",
          "tidyterra","scales","paletteer","plotbiomes",
          "viridis")
install.packages(pkgs)
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE)


output_directory<-"~/Dropbox/bgc_work/stem_bgc/model_school/"
setwd(output_directory)

header<-c("ws.soilw", "wf.canopyw_evap",
 "wf.soilw_outflow", "cs.leafc",
 "cs.frootc", "cs.livestemc",
 "cs.deadstemc", "cs.livecrootc",
 "cs.deadcrootc", "epv.proj_lai",
 "epv.gl_s_sun", "epv.gl_s_shade",
 "epv.m_Kl", "psn_sun.g",
 "psn_sun.Ci", "psn_sun.Ca",
 "psn_sun.A", "psn_shade.g",
 "psn_shade.Ci", "psn_shade.Ca",
 "psn_shade.A")
#data import and basic variables

sim_output_dailiy<-read_tsv("IG_600_21_A_ebf.dayout.ascii",col_names = header)

ggplot(sim_output_dailiy)+geom_point(aes(x=ws.soilw,y=epv.gl_s_sun))

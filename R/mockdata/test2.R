# test for caterogical data

library(here)
library(tidyverse)
library(ape)


dat <- read.csv(here("R","mockdata","mock_data.csv"))

tree <- read.nexus(here("R","mockdata","mock_data_tree.nex"))



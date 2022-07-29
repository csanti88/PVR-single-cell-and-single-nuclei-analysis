#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 14:54:22 2021

@author: megan
"""

#combine loom files
pip install -U loompy
import loompy
files_sn = ["/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_IIT4C.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_9NDRR.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_09ZMR.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_X6H79.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_JV059.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_WSGEO.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_8QQP9.loom"]

files_sc = ["/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_QIGQR.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_LPJH2.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_638JL.loom"]

files_comb = ["/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_IIT4C.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_9NDRR.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_09ZMR.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_X6H79.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_JV059.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_WSGEO.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_8QQP9.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_QIGQR.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_LPJH2.loom",
         "/Users/megan/Documents/Blackshaw_Lab/R/PVR/Raw_Data/BAM/possorted_genome_bam_638JL.loom"]

output_sn_filename = "/Volumes/MYG_5TB/Pvr/analysis/sc_v_sn_manuscript/RNAVelocity/sn.loom"
output_sc_filename = "/Volumes/MYG_5TB/Pvr/analysis/sc_v_sn_manuscript/RNAVelocity/sc.loom"
output_comb_filename = "/Volumes/MYG_5TB/Pvr/analysis/sc_v_sn_manuscript/RNAVelocity/comb.loom"

loompy.combine(files_sn, output_sn_filename, key="Accession")
loompy.combine(files_sc, output_sc_filename, key="Accession")
loompy.combine(files_comb, output_comb_filename, key="Accession")


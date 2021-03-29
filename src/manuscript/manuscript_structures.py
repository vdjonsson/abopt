import pandas as pd
import os 
import antibody_pipeline as ap
import structure 
import energy     


''' This file will run all the structural calculations for the manuscript '''


''' Each individual antibody should be run in parallel in the order below  '''

''' Repair antibody/virus structures, with and without virus and antibody structures removed ''' 
cmd = 'repair_antibodies.py'

''' Run mutational scanning on virus and calculate binding energy from this '''
cmd = 'virus_scan.py'

''' Then merge all mutational scanning data ''' 
cmd = 'merge.py'





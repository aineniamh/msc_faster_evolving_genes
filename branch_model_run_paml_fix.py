import os
from Bio.Phylo.PAML import codeml
            
#indir = 'msa_files/new_dup/'
#g = os.walk(indir)
# = 0
#for r, d, f in g:
    #if os.path.exists(r+'/'+'muscle_alignment.fa'):
    #if r == 'duplicability/msa_files/new_dup/ENSGT00390000002253_ENSMMUP00000016573_ENSNLEP00000018273':
        #print(r)
        #tree_files = []
        #for my_f in f:
            #if my_f.startswith('tree_file.seq_1.'):
                #tree_files.append(my_f)
        #print(tree_files)
        #for t in tree_files:
            #if not os.path.exists(r + '/branch_codeml_results.'+str(t)+'.out'):
            #print(r)
t = 'tree_file.seq_1.nw'
cml = codeml.Codeml()
cml.alignment = 'phylip_alignment.phy'
cml.out_file = 'branch_codeml_results.'+str(t)+'.out'
cml.tree = t

cml.read_ctl_file('/Users/otoolea2/paml4.8/codeml.ctl')

                #cml.print_options()

res = cml.run(command = '/Users/otoolea2/paml4.8/codeml', verbose = True, parse = True)

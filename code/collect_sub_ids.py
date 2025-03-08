import argparse

parser = argparse.ArgumentParser(description=r'finds subsystem proteins to be added to the clustered database')
parser.add_argument('-r', '--rep', action='store', type=str, help='IDs of representative sequences')
parser.add_argument('-s', '--sbs', action='store', type=str, help='subsystem file') 
parser.add_argument('-a', '--allid', action='store', type=str, help='file with all reference IDs to test for non-existing sequences')
parser.add_argument('-o', '--out', action='store', type=str, help='output file with IDs of missed subsystem proteins')
args = parser.parse_args()

outfile = args.out
subfile = args.sbs
repfile =  args.rep
allfile = args.allid

sub_ids = list() # list of IDs in the subsystem table
rep_list = list() # list of IDs in the representative file
all_list = list() # list of IDs in the concatenated FAA file

print(f'Read file with all sequences: {allfile}')
with open(allfile, 'r') as allid:
    for record in allid.readlines():
        all_list.append( record.rstrip() )

print(f'Read file with representative sequences: {repfile}')
with open(repfile, 'r') as rep:
    for record in rep.readlines():
        rep_list.append( record.rstrip() )

print(f'Read subsystem table {subfile}')
with open(subfile, 'r') as sub:
    for line in sub.readlines():
        line = line.split('\t')
        sub_ids.append( line[0] )

print('Calculate difference between subsystems and representatives')
diff = list(set(sub_ids) - set(rep_list))
not_in_all = list(set(diff) - set(all_list))
add = list(set(diff) - set(not_in_all))

with open('mmseqs/not_in_all.txt', 'w') as nia:
    nia.write('\n'.join(set(sub_ids) - set(all_list)))

with open(outfile, 'w') as out:
    out.write('\n'.join(add))
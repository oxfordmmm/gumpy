import gumpy

reference=gumpy.Genome(genbank_file='/Users/fowler/packages/piezo/config/H37rV_v3.gbk',show_progress_bar=True,name="H37rV_v3")

print(reference)

for i in range(1,10):
    reference.save_pickle("H37rV_v3_"+str(i)+".pkl",compression=True,compresslevel=i)

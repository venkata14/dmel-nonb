import os
import concurrent.futures
from time import sleep

rGquadFolder = "R_gquad"
indContigsFolder = "ind-contigs"

if not os.path.exists(rGquadFolder):
    os.makedirs(rGquadFolder)

for n in os.listdir(indContigsFolder):
    if '.fasta' not in n:
        print('Unknown file:', n)
        continue

    fileContigName = n[:-6]
    csvFile = ''.join([fileContigName, '.csv'])
    rFile = ''.join([fileContigName, '.R'])

    fastaFile = os.path.join(indContigsFolder, n)
    fullCsvFile = os.path.join(rGquadFolder, csvFile)
    fullRFile = os.path.join(rGquadFolder, rFile)

    rFileCommand = """library(gquad)


df_gquad <- gquad(x='{0}', xformat='fasta')
write.csv(df_gquad, "gquad_{1}")

df_aphased <- aphased(x='{0}', xformat='fasta')
write.csv(df_aphased, "aphased_{1}")

df_hdna <- hdna(x='{0}', xformat='fasta')
write.csv(df_hdna, "hdna_{1}")

df_slipped <- slipped(x='{0}', xformat='fasta')
write.csv(df_slipped, "slipped_{1}")

df_str <- str(x='{0}', xformat='fasta')
write.csv(df_str, "str_{1}")

df_tfo <- tfo(x='{0}', xformat='fasta')
write.csv(df_tfo, "tfo_{1}")

df_zdna <- zdna(x='{0}', xformat='fasta')
write.csv(df_zdna, "zdna_{1}")


print('Done With: {2}')""".format(fastaFile, csvFile, fullRFile) # fullCsvFile is not needed as this R file is in the same directory as the csv file

    with open(fullRFile, 'w') as r:
        r.write(rFileCommand)
    

def runRFile(command):
    os.system("Rscript {}".format(command))

rFiles = []

for n in os.listdir(rGquadFolder):
    if '.R' not in n:
        continue
    rFiles.append(os.path.join(rGquadFolder, n))

executor = concurrent.futures.ProcessPoolExecutor(20)
futures = [executor.submit(runRFile, file) for file in rFiles]
concurrent.futures.wait(futures)

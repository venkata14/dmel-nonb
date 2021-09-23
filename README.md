# dmel-nonb

This repository has all the code for the paper ***

Extraction of IPD Values for *D. melanogaster*
---

1. Initialize the environment 

```
sudo apt-get update
sudo apt-get install build-essential wget locales
sudo dpkg-reconfigure locales
sudo apt-get install rsync
sudo locale-gen en_US.UTF-8
```

The tool used for the alignment and extraction was `smrttools 7.0.0`

The download can be found here: https://www.pacb.com/support/software-downloads/
The download instructions can be found here: https://www.pacb.com/wp-content/uploads/SMRT_Link_Installation_v701.pdf

2. Install the datasets. Datasets from (ref)

Install the datasets
```
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro1_24NOV2013_398.tgz
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro2_25NOV2013_399.tgz
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro3_26NOV2013_400.tgz
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro4_28NOV2013_401.tgz
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro5_29NOV2013_402.tgz
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro6_1DEC2013_403.tgz
```

Getting md5sums

```
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro1_24NOV2013_398.md5
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro2_25NOV2013_399.md5
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro3_26NOV2013_400.md5
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro4_28NOV2013_401.md5
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro5_29NOV2013_402.md5
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro6_1DEC2013_403.md5
```

Verify archives have been downloaded properly
```
md5sum -c Dro1_24NOV2013_398.md5 > Dro1_24NOV2013_398.md5.checkresult
md5sum -c Dro2_25NOV2013_399.md5 > Dro2_25NOV2013_399.md5.checkresult
md5sum -c Dro3_26NOV2013_400.md5 > Dro3_26NOV2013_400.md5.checkresult
md5sum -c Dro4_28NOV2013_401.md5 > Dro4_28NOV2013_401.md5.checkresult
md5sum -c Dro5_29NOV2013_402.md5 > Dro5_29NOV2013_402.md5.checkresult
md5sum -c Dro6_1DEC2013_403.md5  > Dro6_1DEC2013_403.md5.checkresult
```

3. Extract files with `tar -xvzf <file.tgz>`

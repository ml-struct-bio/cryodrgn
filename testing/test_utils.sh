python ../cryodrgn/commands/downsample.py data/hand.mrcs -o output/hand.mrcs -D 32 
python ../cryodrgn/commands/downsample.py data/toymodel_small_nocenter.mrc -D 24 --is-vol -o output/test.mrc
python ../cryodrgn/commands/downsample.py data/hand.mrcs -o output/hand.mrcs --chunk 10 -D 32

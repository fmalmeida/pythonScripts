# make dir
mkdir -p build

# build conda
conda build --output-folder build/ conda.recipe

# convert to osx
conda convert -p osx-64 $(find build -name "falmeida-py*.tar.bz2")

# upload osx
anaconda upload $(find osx-64 -name "falmeida-py*.tar.bz2") --force

# rm dirs
rm -rf build osx-64

# save new help message
bash build_helps.sh

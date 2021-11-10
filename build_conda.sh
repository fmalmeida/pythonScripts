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
python3 falmeida-py-runner.py -h &> docs/help_message.txt
python3 falmeida-py-runner.py tsv2markdown -h &> docs/tsv2markdown_help.txt
python3 falmeida-py-runner.py splitgbk -h &> docs/splitgbk_help.txt
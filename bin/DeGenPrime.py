import argparse

#creating an ArgParse to accept file input
parser = argparse.ArgumentParser()
required = parser.add_argument_group("Required Arguments")
required.add_argument('-i','--input',type = str, help='path to intended input file')
required.add_argument('-o','--output',type = str, help='path to intended output file')
args = parser.parse_args()
file_in = open(args.input)


#program will split here - need to check to see if the file is aligned or not
# if file is aligned, continue on with analysis
# if file is not aligned, send to widget


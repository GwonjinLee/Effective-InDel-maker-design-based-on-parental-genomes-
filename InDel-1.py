import sys

inFile = open(sys.argv[1], "r")
outFile = open(sys.argv[2], "w")

for line in inFile:
    line = line.strip("\n")  # remove new lines

    if line.startswith('#'):
        continue
    else:
        List = line.split()  # Split the line into a list
        length = len(List[4]) - len(List[3])

        if 85 > length > 25 or -85 < length < -25:
            length = len(List[4]) - len(List[3])
            outFile.write('{}\t{}\n'.format(line, length))

# Close the files
inFile.close()
outFile.close()
